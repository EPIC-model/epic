!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : pi, f12, max_num_parcels
    use parcel_container, only : parcels, n_parcels, get_delx
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, nz
    use options, only : parcel

    implicit none

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:),kc2(:)
    integer :: loca(max_num_parcels)
    integer :: node(max_num_parcels)

    !Other variables:
    double precision:: vmin, delx,delz,dsq,dsqmin,x_small,z_small
    integer:: i,ic,i0, ib,k,m,j, is
    integer:: ix,iz,ix0,iz0

    public :: find_nearest

    contains

        ! if ibig(n) is zero, no parcel found for isma(n)
        subroutine find_nearest(isma, ibig, nmerge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: ibig(:)
            integer, intent(out) :: nmerge

            if (.not. allocated(nppc)) then
                allocate(nppc(ncell))
                allocate(kc1(ncell))
                allocate(kc2(ncell))
            endif

            vmin = vcell / dble(parcel%vmin_fraction)

            nmerge=0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc=0 !nppc(ic) will contain the number of parcels in grid cell ic

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do i = 1, n_parcels
                ix=int(dxi(1)*(parcels%position(i,1)-lower(1)))
                iz=int(dxi(2)*(parcels%position(i,2)-lower(2)))

                ! Cell index of parcel:
                ic=1+ix+nx*iz !This runs from 1 to ncell

                ! Accumulate number of parcels in this grid cell:
                nppc(ic)=nppc(ic)+1

                ! Store grid cell that this parcel is in:
                loca(i)=ic

                if (parcels%volume(i) < vmin) then
                    nmerge = nmerge + 1
                endif
            enddo

            if (nmerge == 0) then
                return
            endif

            ! allocate arrays
            allocate(isma(0:nmerge))
            allocate(ibig(nmerge))

            isma = 0
            ibig = 0
            nmerge = 0

            ! Find arrays kc1(ic) & kc2(ic) which indicate the parcels in grid cell ic
            ! through i = node(k), for k = kc1(ic),kc2(ic):
            kc1(1)=1
            do ic=1,ncell-1
                kc1(ic+1)=kc1(ic)+nppc(ic)
            enddo

            kc2=kc1-1
            do i=1,n_parcels
                ic=loca(i)
                k=kc2(ic)+1
                node(k)=i
                kc2(ic)=k

                if (parcels%volume(i) < vmin) then
                    nmerge = nmerge + 1
                    isma(nmerge) = i
                endif
                ! reset loca (for use later)
                loca(i) = 0
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:
            j = 0
            ! j counts the actual number of mergers found
            do m=1,nmerge
                i0=isma(m)
                x_small=parcels%position(i0,1)
                z_small=parcels%position(i0,2)
                ! Parcel i0 is small and should be merged; find closest other:
                ix0=mod(nint(dxi(1)*(x_small-lower(1))),nx) ! ranges from 0 to nx-1
                iz0=nint(dxi(2)*(z_small-lower(2)))         ! ranges from 0 to nz

                ! Grid point (ix0,iz0) is closest to parcel i0

                dsqmin = product(extent)
                ib = 0

!                 print *, max(0,iz0-1), min(nz,iz0+1), z_small, lower(2), dxi(2)

                ! Loop over 8 cells surrounding (ix0,iz0):
                do iz=max(0,iz0-1),min(nz-1,iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do ix=ix0-1,ix0
                        ! Cell index (accounting for x periodicity):
                        ic=1+mod(nx+ix,nx)+nx*iz
                        ! Search parcels for closest:
                        do k=kc1(ic),kc2(ic)
                            i = node(k)
                            if (parcels%volume(i) >= parcels%volume(i0)) then
                                if (i .ne. i0) then
                                    delz = parcels%position(i,2) - z_small
                                    delx = get_delx(parcels%position(i,1), x_small) ! works across periodic edge
                                    ! Minimise dsqmin
                                    dsq= delz * delz + delx * delx
                                    if (dsq < dsqmin) then
        !                                 print *, "i = ", i
                                        dsqmin = dsq
                                        ib = i
                                    endif
                                endif
                            endif
                        enddo
                    enddo
                enddo
                if (ib > 0) then
                    ! Store the index of the parcel to be merged with:
                    j = j + 1
                    isma(j) = i0
                    ibig(j) = ib
                    loca(i0) = loca(i0) + 1
!                     print *, "ib = ", ib
                    loca(ib) = loca(ib) + 1
                endif
            enddo
            ! Actual total number of mergers:
            nmerge = j

            j = 0
            do m = 1, nmerge
                ib = ibig(m)
                is = isma(m)
                if ((loca(ib) > 1) .and. (loca(is) > 1)) then
                    ! remove link between "is" and "ib"
                    loca(ib) = loca(ib) - 1
                    loca(is) = loca(is) - 1
!                     print *, "remove link"
!                     stop
                else
                    j = j + 1
                    isma(j) = isma(m)
                    ibig(j) = ibig(m)
                endif
            enddo
            nmerge = j


            do m = 1, nmerge
                ib = ibig(m)
                is = isma(m)
                if (loca(ib) < loca(is)) then
                    isma(m) = ib
                    ibig(m) = is
                endif
            enddo

!             print *, isma(1:nmerge)
!             print *, ibig(1:nmerge)
!             stop

        end subroutine find_nearest

end module parcel_nearest
