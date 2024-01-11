module parcel_mixing
    use constants, only : pi, f12
    use parcel_ops, only : get_delx
    use parcel_container, only : parcels, n_parcels
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_bot_parcels, bot_parcels     &
                                       , n_top_parcels, top_parcels
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, nz, vmin, max_num_parcels, lmin
    use options, only : parcel
    use timer, only : start_timer, stop_timer

    implicit none

    integer:: mixing_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: top_nppc(:), bot_nppc(:)!, kc1(:),kc2(:)
    integer, allocatable :: top_loca(:), bot_loca(:)
    integer, allocatable :: top_pid(:), bot_pid(:)

!     integer, allocatable :: node(:)

    !Other variables:
!     double precision:: delx,delz,dsq,dsqmin,x_small,z_small
!     integer:: ic,is,i,k,m,j, n
!     integer:: ix,ix0,iz0

    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels

            ! 1. find small parcels in cells 0 and nz
            ! 2. find closest parcels
            ! 3. mix properties

            call find_interior_parcels

            call surface_to_interior(0,    n_bot_parcels, bot_parcels, top_nppc, top_loca, top_pid)
            call surface_to_interior(nz+1, n_top_parcels, top_parcels, bot_nppc, bot_loca, bot_pid)

        end subroutine mix_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine surface_to_interior(iz, n_spar, spar, nppc, loca, pid)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_spar
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: nppc(:)
            integer,                             intent(in)    :: loca(:)
            integer,                             intent(in)    :: pid(:)
            integer                                            :: n_interior, n, ix, i, m, ix0, ic, k
            double precision                                   :: length, xs, zs, delx, delz
            integer, allocatable                               :: isma(:)
            integer, allocatable                               :: iclo(:)

            call start_timer(mixing_timer)

            m = 1
            do n = 1, n_spar
                length = get_surface_parcel_length(n, spar)
                if (length < lmin) then

                    ix = mod(int(dxi(1) * (spar%position(n) - lower(1))), nx)

                    ! Cell index of parcel:
                    i = 1 + ix !This runs from 1 to nx

                    ! Accumulate number of parcels in this grid cell:
                    nppc(i) = nppc(i) + 1

                    ! Store grid cell that this parcel is in:
                    loca(m) = i

                    ! Store parcel index:
                    pid(m) = n

                    m = m + 1
                endif
            enddo

            ! number of small surface parcels
            nmix = m - 1


            ! -----------------------------------------------------------------
            ! Find closest interior parcel to small surface parcel:

            allocate(isma(nmix))
            allocate(iclo(nmix))

            isma = 0
            iclo = 0

            do m = 1, nmix

                ! grid cell this parcel is in:
                ix0 = loca(m)

                ! parcel index
                is = pid(m)

                ! position of surface parcel
                xs = spar%position(is)
                zs = lower(2) + dx(2) * dble(iz)

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





!             call find_parcels
!
!             if (nmerge == 0) then
!                 call stop_timer(merge_nearest_timer)
!                 return
!             endif
!
!             ! allocate arrays
!             allocate(isma(0:nmerge))
!             allocate(iclo(nmerge))
!
!             isma = 0
!             iclo = 0
!
!             ! Find arrays kc1(ij) & kc2(ij) which indicate the parcels in grid cell ij
!             ! through n = node(k), for k = kc1(ij),kc2(ij):
!             kc1(1) = 1
!             do ij = 1, ncell-1
!                 kc1(ij+1) = kc1(ij) + nppc(ij)
!             enddo
!
!             kc2 = kc1 - 1
!             j = 0
!             do n = 1, n_parcels
!                 ij = loca(n)
!                 k = kc2(ij) + 1
!                 node(k) = n
!                 kc2(ij) = k
!
!                 if (parcels%volume(n) < vmin) then
!                     j = j + 1
!                     isma(j) = n
!                 endif
!             enddo
!
!             !---------------------------------------------------------------------
!             ! Now find the nearest grid point to each small parcel (to be merged)
!             ! and search over the surrounding 8 grid cells for the closest parcel:
!
!             ! SB: do not use temporary (j) index here, so we will be able to parallelise.
!             ! Rather, stop if no nearest parcel found  in surrounding grid boxes
!             do m = 1, nmerge
!                 is = isma(m)
!                 x_small = parcels%position(1, is)
!                 z_small = parcels%position(2, is)
!                 ! Parcel "is" is small and should be merged; find closest other:
!                 ix0 = mod(nint(dxi(1) * (x_small - lower(1))), nx) ! ranges from 0 to nx-1
!                 iz0 = nint(dxi(2) * (z_small - lower(2)))          ! ranges from 0 to nz
!
!                 ! Grid point (ix0,iz0) is closest to parcel "is"
!
!                 dsqmin = product(extent)
!                 ic = 0
!
!                 ! Loop over 8 cells surrounding (ix0,iz0):
!                 do iz = max(0, iz0-1), min(nz-1, iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
!                     do ix = ix0-1, ix0
!                         ! Cell index (accounting for x periodicity):
!                         ij = 1 + mod(nx + ix, nx) + nx * iz
!                         ! Search small parcels for closest other:
!                         do k = kc1(ij), kc2(ij)
!                             n = node(k)
!                             if (n .ne. is) then
!                                 delz = parcels%position(2, n) - z_small
!                                 if (delz*delz < dsqmin) then
!                                     delx = get_delx(parcels%position(1, n), x_small) ! works across periodic edge
!                                     ! Minimise dsqmin
!                                     dsq = delz * delz + delx * delx
!                                     if (dsq < dsqmin) then
!                                         dsqmin = dsq
!                                         ic = n
!                                     endif
!                                 endif
!                             endif
!                         enddo
!                     enddo
!                 enddo
!
!                 if (ic==0) then
!                   print *, 'Merge error: no near neighbour found.'
!                   stop
!                 end if
!
!                 ! Store the index of the parcel to be potentially merged with:
!                 isma(m) = is
!                 iclo(m) = ic
!             enddo

            call stop_timer(mixing_timer)

        end subroutine surface_to_interior

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine find_interior_parcels
            integer :: k, i, n, ix, ib, it

            if (.not. allocated(top_nppc)) then
                allocate(top_nppc(nx))
                allocate(bot_nppc(nx))
!                 allocate(kc1(nx))
!                 allocate(kc2(nx))
                allocate(top_loca(max_num_parcels / nz))
                allocate(bot_loca(max_num_parcels / nz))

                allocate(top_pid(max_num_parcels / nz))
                allocate(bot_pid(max_num_parcels / nz))

!                 allocate(node(max_num_parcels))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            top_nppc = 0 !top_nppc(i) will contain the number of parcels in grid cell i
            bot_nppc = 0

            ib = 1
            it = 1

            ! Bin interior parcels in cells:
            do n = 1, n_parcels

                ! we only consider big parcels
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                k = nint(dxi(2) * (parcels%position(2, n) - lower(2)))

                ix = mod(int(dxi(1) * (parcels%position(1, n) - lower(1))), nx)

                ! Cell index of parcel:
                i = 1 + ix !This runs from 1 to nx

                if (k == nz) then
                    !
                    ! top surface
                    !

                    ! Accumulate number of parcels in this grid cell:
                    top_nppc(i) = top_nppc(i) + 1

                    ! Store grid cell that this parcel is in:
                    top_loca(it) = i

                    ! Store parcel index:
                    top_pid(it) = n

                    it = it + 1

                else if (k == 0) then
                    !
                    ! bottom surface
                    !
                    bot_nppc(i) = bot_nppc(i) + 1

                    bot_loca(ib) = i

                    bot_pid(ib) = n

                    ib = ib + 1

                endif
            enddo

        end subroutine find_interior_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mixing(spar, isma, iclo, nmix)
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: nmix
            integer                                            :: l, m, n, ic, is
            double precision                                   :: ai, vmix

            loca = zero

            !------------------------------------------------------------------
            ! Figure out weights:

            l = 0
            do m = 1, nmix

                ! Index of closest interior parcel
                ic = iclo(m)

                if (loca(ic) == 0) then
                    ! Start a new mixing parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    wm(l) = parcels%volume(ic) / vcell

                    buoym(l) = wm(l) * parcels%buoyancy(ic)

                    vortm(l) = wm(l) * parcels%vorticity(ic)

                endif

                ! Sum up all the small surface parcel contributions:
                is = isma(m)
                n = loca(ic)

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
            loca = zero
            l = 0

            do m = 1, nmix
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    parcels%buoyancy(ic) = buoym(l)
                    parcels%vorticity(ic) = vortm(l)
                endif

                is = isma(m)
                n = loca(ic)

                spar%buoyancy(is) = buoym(l)
                spar%vorticity(is) = vortm(l)
            enddo

        end subroutine mixing

end module parcel_mixing
