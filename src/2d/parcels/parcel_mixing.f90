module parcel_mixing
    use constants, only : one, zero, two
    use parcel_ops, only : get_delx
    use parcel_container, only : parcels, n_parcels
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_bot_parcels, bot_parcels     &
                                       , n_top_parcels, top_parcels     &
                                       , get_surface_parcel_length
    use parameters, only : dx, dxi, vcell, hli, lower, extent, lcell    &
                         , nx, nz, vmin, max_num_parcels, lmin          &
                         , max_num_surf_parcels
    use options, only : parcel
    use timer, only : start_timer, stop_timer

    implicit none

    integer :: mixing_timer

    private

    !Used for searching for possible parcel mixer:
    integer, allocatable :: nppc(:), kc1(:), kc2(:)
    integer, allocatable :: loca(:)
    integer, allocatable :: node(:)
    logical, allocatable :: top_mixed(:), bot_mixed(:)

    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels

            call start_timer(mixing_timer)

            ! 1. find small parcels in cells 0 and nz
            ! 2. find closest parcels
            ! 3. mix properties

            if (.not. allocated(top_mixed)) then
                allocate(top_mixed(max_num_surf_parcels))
                allocate(bot_mixed(max_num_surf_parcels))
            endif

            call interior2surface(0,  n_bot_parcels, bot_parcels, bot_mixed)
            call interior2surface(nz, n_top_parcels, top_parcels, top_mixed)

            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
            endif

            call find_interior_parcels

            call surface2interior(0,  n_bot_parcels, bot_parcels, bot_mixed)
            call surface2interior(nz, n_top_parcels, top_parcels, top_mixed)

            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
            endif


            call stop_timer(mixing_timer)

        end subroutine mix_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine surface2interior(iz, n_spar, spar, is_mixed)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_spar
            type(surface_parcel_container_type), intent(inout) :: spar
            logical, allocatable,                intent(in)    :: is_mixed(:)
            integer                                            :: n, ix, ij, m, ix0, ic, k, is
            double precision                                   :: xs, zs, delx, delz, dsq, dsqmin
            integer, allocatable                               :: isma(:)
            integer, allocatable                               :: iclo(:)
            integer                                            :: nmix

            ! -----------------------------------------------------------------
            ! Number of small surface parcels:
            nmix = 0
            do n = 1, n_spar
                if ((spar%volume(n) < vmin) .and. (.not. is_mixed(n))) then
                    nmix = nmix + 1
                endif
            enddo

            if (nmix == 0) then
                return
            endif

            ! -----------------------------------------------------------------
            ! Setup search arrays:

            allocate(isma(nmix))
            allocate(iclo(nmix))

            isma = 0
            iclo = 0

            ! Fill isma array:
            m = 0
            do n = 1, n_spar
                if ((spar%volume(n) < vmin) .and. (.not. is_mixed(n))) then
                    m = m + 1
                    isma(m) = n
                endif
            enddo

            ! -----------------------------------------------------------------
            ! Find closest interior parcel to small surface parcel:

            do m = 1, nmix

                is = isma(m)

                ! position of surface parcel
                xs = spar%position(is)
                zs = lower(2) + dx(2) * dble(iz)

                ! grid cell this parcel is in:
                ix0 = mod(nint(dxi(1) * (xs - lower(1))), nx)

                dsqmin = two * dx(1) * dx(2)
                ic = 0

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do ix = ix0-1, ix0
                    ! Cell index (accounting for x periodicity):
                    ij = 1 + mod(nx + ix, nx) + nx * min(iz, 1)
                    ! Search small parcels for closest other:
                    do k = kc1(ij), kc2(ij)
                        n = node(k)
                        delz = parcels%position(2, n) - zs
                        if (delz * delz < dsqmin) then
                            delx = get_delx(parcels%position(1, n), xs) ! works across periodic edge
                            ! Minimise dsqmin
                            dsq = delz * delz + delx * delx
                            if (dsq < dsqmin) then
                                dsqmin = dsq
                                ic = n
                            endif
                        endif
                    enddo
                enddo

                ! Store the index of the parcel to be mixed with:
                isma(m) = is
                iclo(m) = ic

                ! If ic == 0, then no near interior parcel is found. We must ignore this mix.
                if (ic == 0) then
                    isma(m) = 0
                    nmix = nmix - 1
                endif
            enddo

            isma = pack(isma, isma /= 0)
            iclo = pack(iclo, iclo /= 0)


            ! -----------------------------------------------------------------
            ! Mix small surface parcel with closest interior parcel:

            call mix_surface2interior(spar, isma, iclo, nmix)

            ! -----------------------------------------------------------------
            ! Final clean-up:

            deallocate(isma)
            deallocate(iclo)

        end subroutine surface2interior

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine find_interior_parcels
            integer :: iz, ij, n, ix, k

            if (.not. allocated(nppc)) then
                allocate(nppc(2 * nx))
                allocate(kc1(2 * nx))
                allocate(kc2(2 * nx))
                allocate(loca(max_num_parcels))
                allocate(node(max_num_parcels))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(i) will contain the number of parcels in grid cell i

            ! Bin interior parcels in cells:
            do n = 1, n_parcels

                ! we only consider big parcels
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                iz = nint(dxi(2) * (parcels%position(2, n) - lower(2)))

                if ((iz == nz) .or. (iz == 0)) then
                    !
                    ! top or bottom surface
                    !

                    ix = mod(int(dxi(1) * (parcels%position(1, n) - lower(1))), nx)

                    ! Cell index of parcel:
                    ij = 1 + mod(nx + ix, nx) + nx * min(iz, 1)


                    ! Accumulate number of parcels in this grid cell:
                    nppc(ij) = nppc(ij) + 1

                    ! Store grid cell that this parcel is in:
                    loca(n) = ij
                endif
            enddo

            ! Find arrays kc1(i) & kc2(i) which indicate the parcels in grid cell i
            ! through n = node(k), for k = kc1(i),kc2(i):
            kc1(1) = 1
            do ij = 1, 2 * nx - 1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            do n = 1, n_parcels
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                iz = nint(dxi(2) * (parcels%position(2, n) - lower(2)))
                if ((iz == nz) .or. (iz == 0)) then
                    ij = loca(n)
                    k = kc2(ij) + 1
                    node(k) = n
                    kc2(ij) = k
                endif
            enddo

        end subroutine find_interior_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mix_surface2interior(spar, isma, iclo, nmix)
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: nmix
            integer                                            :: l, m, n, ic, is
            double precision                                   :: buoym(nmix), vortm(nmix), vm(nmix), vmix
            integer                                            :: lclo(n_parcels)

            lclo = zero


            !------------------------------------------------------------------
            ! Figure out weights:

            l = 0
            do m = 1, nmix

                ! Index of closest interior parcel
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    ! Start a new mixing parcel, indexed l:
                    l = l + 1
                    lclo(ic) = l

                    vm(l) = parcels%volume(ic)


                    buoym(l) = vm(l) * parcels%buoyancy(ic)

                    vortm(l) = vm(l) * parcels%vorticity(ic)

                endif

                ! Sum up all the small surface parcel contributions:
                is = isma(m)
                n = lclo(ic)

                vm(n) = vm(n) + spar%volume(is)

                buoym(n) = buoym(n) + spar%volume(is) * spar%buoyancy(is)

                vortm(n) = vortm(n) + spar%volume(is) * spar%vorticity(is)
            enddo

            !------------------------------------------------------------------
            ! Obtain the mixed values:
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmix = one / vm(m)

                buoym(m) = vmix * buoym(m)

                vortm(m) = vmix * vortm(m)
            enddo


            !------------------------------------------------------------------
            ! Apply the blended values:
            lclo = zero
            l = 0

            do m = 1, nmix
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    l = l + 1
                    lclo(ic) = l

                    parcels%buoyancy(ic) = buoym(l)
                    parcels%vorticity(ic) = vortm(l)
                endif

                is = isma(m)
                n = lclo(ic)

                spar%buoyancy(is) = buoym(l)
                spar%vorticity(is) = vortm(l)
            enddo

        end subroutine mix_surface2interior

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine interior2surface(iz, n_spar, spar, is_mixed)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_spar
            type(surface_parcel_container_type), intent(inout) :: spar
            logical, allocatable,                intent(inout) :: is_mixed(:)
            integer                                            :: n, ix, i, m, ix0, ic, k, is
            double precision                                   :: xs, zs, delx, delz, dsq, dsqmin
            integer, allocatable                               :: isma(:)
            integer, allocatable                               :: iclo(:)
            integer                                            :: nmix, ii, ij, l

            is_mixed = .false.

            ! -----------------------------------------------------------------
            ! Number of small interior parcels:
            nmix = 0
            l = 0
            do n = 1, n_parcels
                if (parcels%volume(n) < vmin) then
                    l = l + 1
                    ii = nint(dxi(2) * (parcels%position(2, n) - lower(2)))
                    if (ii == iz) then
                        nmix = nmix + 1
                    endif
                endif
            enddo

            if (nmix == 0) then
                return
            endif

            if (.not. allocated(nppc)) then
                allocate(nppc(nx))
                allocate(kc1(nx))
                allocate(kc2(nx))
                allocate(loca(max_num_surf_parcels))
                allocate(node(max_num_surf_parcels))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(i) will contain the number of parcels in grid cell i

            ! Bin interior parcels in cells:
            do n = 1, n_spar

                ix = mod(int(dxi(1) * (spar%position(n) - lower(1))), nx)

                ! Cell index of parcel:
                ij = 1 + ix

                ! Accumulate number of parcels in this grid cell:
                nppc(ij) = nppc(ij) + 1

                ! Store grid cell that this parcel is in:
                loca(n) = ij
            enddo

            ! Find arrays kc1(i) & kc2(i) which indicate the parcels in grid cell i
            ! through n = node(k), for k = kc1(i),kc2(i):
            kc1(1) = 1
            do ij = 1, nx - 1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            do n = 1, n_spar
                ij = loca(n)
                k = kc2(ij) + 1
                node(k) = n
                kc2(ij) = k
            enddo

            ! -----------------------------------------------------------------
            ! Setup search arrays:

            allocate(isma(nmix))
            allocate(iclo(nmix))

            isma = 0
            iclo = 0

            ! Fill isma array:
            m = 0
            do n = 1, n_parcels
                if (parcels%volume(n) < vmin) then
                    ii = nint(dxi(2) * (parcels%position(2, n) - lower(2)))
                    if (ii == iz) then
                        m = m + 1
                        isma(m) = n
                    endif
                endif
            enddo

            ! -----------------------------------------------------------------
            ! Find closest interior parcel to small surface parcel:

            do m = 1, nmix

                is = isma(m)

                ! position of surface parcel
                xs = parcels%position(1, is)
                zs = parcels%position(2, is)

                ! grid cell this parcel is in:
                ix0 = mod(nint(dxi(1) * (xs - lower(1))), nx)

                dsqmin = two * dx(1) * dx(2)
                ic = 0

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do ix = ix0-1, ix0
                    ! Cell index (accounting for x periodicity):
                    i = 1 + mod(nx + ix, nx)
                    ! Search small parcels for closest other:
                    do k = kc1(i), kc2(i)
                        n = node(k)
                        delz = lower(2) + dx(2) * dble(iz) - zs
                        if (delz * delz < dsqmin) then
                            delx = get_delx(spar%position(n), xs) ! works across periodic edge
                            ! Minimise dsqmin
                            dsq = delz * delz + delx * delx
                            if (dsq < dsqmin) then
                                dsqmin = dsq
                                ic = n
                            endif
                        endif
                    enddo
                enddo

                if (ic == 0) then
                    print *, 'Mixing error: no near surface parcel found.'
                    stop
                endif

                is_mixed(ic) = .true.

                ! Store the index of the parcel to be mixed with:
                isma(m) = is
                iclo(m) = ic
            enddo


            ! -----------------------------------------------------------------
            ! Mix small interor parcel with closest surface parcel:

            call mix_interior2surface(spar, isma, iclo, nmix, n_spar)

            ! -----------------------------------------------------------------
            ! Final clean-up:

            deallocate(isma)
            deallocate(iclo)

        end subroutine interior2surface

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mix_interior2surface(spar, isma, iclo, nmix, n_spar)
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: nmix
            integer,                             intent(in)    :: n_spar
            integer                                            :: l, m, n, ic, is
            double precision                                   :: buoym(nmix), vortm(nmix), vm(nmix), vmix
            integer                                            :: lclo(n_spar)

            lclo = zero

            !------------------------------------------------------------------
            ! Figure out weights:

            l = 0
            do m = 1, nmix

                ! Index of closest interior parcel
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    ! Start a new mixing parcel, indexed l:
                    l = l + 1
                    lclo(ic) = l

                    vm(l) = spar%volume(ic)

                    buoym(l) = vm(l) * spar%buoyancy(ic)

                    vortm(l) = vm(l) * spar%vorticity(ic)

                endif

                ! Sum up all the small surface parcel contributions:
                is = isma(m)
                n = lclo(ic)

                vm(n) = vm(n) + parcels%volume(is)

                buoym(n) = buoym(n) + parcels%volume(is) * parcels%buoyancy(is)

                vortm(n) = vortm(n) + parcels%volume(is) * parcels%vorticity(is)
            enddo

            !------------------------------------------------------------------
            ! Obtain the mixed values:
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmix = one / vm(m)

                buoym(m) = vmix * buoym(m)

                vortm(m) = vmix * vortm(m)
            enddo


            !------------------------------------------------------------------
            ! Apply the blended values:
            lclo = zero
            l = 0

            do m = 1, nmix
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    l = l + 1
                    lclo(ic) = l

                    spar%buoyancy(ic) = buoym(l)
                    spar%vorticity(ic) = vortm(l)
                endif

                is = isma(m)
                n = lclo(ic)

                parcels%buoyancy(is) = buoym(l)
                parcels%vorticity(is) = vortm(l)
            enddo

        end subroutine mix_interior2surface

end module parcel_mixing
