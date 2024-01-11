module parcel_mixing
    use constants, only : one, zero
    use parcel_ops, only : get_delx
    use parcel_container, only : parcels, n_parcels
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_bot_parcels, bot_parcels     &
                                       , n_top_parcels, top_parcels     &
                                       , get_surface_parcel_length
    use parameters, only : dx, dxi, vcell, hli, lower, extent, lcell, nx, nz, vmin, max_num_parcels, lmin
    use options, only : parcel
    use timer, only : start_timer, stop_timer

    implicit none

    integer:: mixing_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:), kc2(:)
    integer, allocatable :: loca(:)
    integer, allocatable :: node(:)

    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels

            ! 1. find small parcels in cells 0 and nz
            ! 2. find closest parcels
            ! 3. mix properties

            call find_interior_parcels

            call surface_to_interior(0,    n_bot_parcels, bot_parcels)
            call surface_to_interior(nz+1, n_top_parcels, top_parcels)

        end subroutine mix_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine surface_to_interior(iz, n_spar, spar)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_spar
            type(surface_parcel_container_type), intent(inout) :: spar
            integer                                            :: n, ix, i, m, ix0, ic, k, is
            double precision                                   :: length, xs, zs, delx, delz, dsq, dsqmin
            integer, allocatable                               :: isma(:)
            integer, allocatable                               :: iclo(:)
            integer                                            :: nmix

            call start_timer(mixing_timer)

            ! -----------------------------------------------------------------
            ! Number of small surface parcels:
            nmix = 0
            do n = 1, n_spar
                length = get_surface_parcel_length(n, spar)
                if (length < lmin) then
                    nmix = nmix + 1
                endif
            enddo

            if (nmix == 0) then
                call stop_timer(mixing_timer)
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
                length = get_surface_parcel_length(n, spar)
                if (length < lmin) then
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
                ix0 = mod(int(dxi(1) * (xs - lower(1))), nx)

                dsqmin = product(extent)
                ic = 0

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do ix = ix0-1, ix0
                    ! Cell index (accounting for x periodicity):
                    i = 1 + mod(nx + ix, nx)
                    ! Search small parcels for closest other:
                    do k = kc1(i), kc2(i)
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

                if (ic == 0) then
                    print *, 'Merge error: no near neighbour found.'
                    stop
                endif

                ! Store the index of the parcel to be mixed with:
                isma(m) = is
                iclo(m) = ic
            enddo


            ! -----------------------------------------------------------------
            ! Mix small surface parcel with closest interior parcel:

            call mixing(spar, isma, iclo, nmix)

            ! -----------------------------------------------------------------
            ! Final clean-up:

            deallocate(isma)
            deallocate(iclo)

            call stop_timer(mixing_timer)

        end subroutine surface_to_interior

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine find_interior_parcels
            integer :: iz, ij, n, ix, l, k

            if (.not. allocated(nppc)) then
                allocate(nppc(2 * nx))
                allocate(kc1(2 * nx))
                allocate(kc2(2 * nx))
                allocate(loca(2 * max_num_parcels / nz))
                allocate(node(2 * max_num_parcels / nz))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(i) will contain the number of parcels in grid cell i

            ! Bin interior parcels in cells:
            l = 1
            do n = 1, n_parcels

                ! we only consider big parcels
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                iz = min(int(dxi(2) * (parcels%position(2, n) - lower(2))), nz-1)

                if ((iz == nz - 1) .or. (iz == 0)) then
                    !
                    ! top or bottom surface
                    !

                    ix = mod(int(dxi(1) * (parcels%position(1, n) - lower(1))), nx)

                    ! Cell index of parcel:
                    ij = 1 + ix + nx * min(iz, 1)


                    ! Accumulate number of parcels in this grid cell:
                    nppc(ij) = nppc(ij) + 1

                    ! Store grid cell that this parcel is in:
                    loca(l) = ij

                    l = l + 1
                endif
            enddo

            ! Find arrays kc1(i) & kc2(i) which indicate the parcels in grid cell i
            ! through n = node(k), for k = kc1(i),kc2(i):
            kc1(1) = 1
            do ij = 1, 2 * nx - 1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            do n = 1, l - 1
                ij = loca(n)
                k = kc2(ij) + 1
                node(k) = n
                kc2(ij) = k
            enddo

        end subroutine find_interior_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mixing(spar, isma, iclo, nmix)
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: nmix
            integer                                            :: l, m, n, ic, is
            double precision                                   :: ai, wmix, length
            double precision                                   :: buoym(nmix), vortm(nmix), wm(nmix)
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

                    wm(l) = parcels%volume(ic) / vcell

                    buoym(l) = wm(l) * parcels%buoyancy(ic)

                    vortm(l) = wm(l) * parcels%vorticity(ic)

                endif

                ! Sum up all the small surface parcel contributions:
                is = isma(m)
                n = lclo(ic)

                length = get_surface_parcel_length(is, spar)

                ai = length / lcell

                wm(n) = wm(n) + ai

                buoym(n) = buoym(n) + ai * spar%buoyancy(is)

                vortm(n) = vortm(n) + ai * spar%vorticity(is)
            enddo

            !------------------------------------------------------------------
            ! Obtain the mixed values:
            do m = 1, l

                wmix = one / wm(m)

                buoym(m) = wmix * buoym(m)

                vortm(m) = wmix * vortm(m)
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

        end subroutine mixing

end module parcel_mixing
