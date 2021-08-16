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
    double precision:: vmin, delx,delz,dsq,dsqmin,x_small,z_small ! ,vmerge,vmergemin
    integer :: i,ic,i0,imin,k,m,j, is, ib
    integer :: ix,iz,ix0,iz0

    public :: find_nearest

    contains

        ! if ibig(n) is zero, no parcel found for isma(n)
        subroutine find_nearest(isma, ibig, nmerge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: ibig(:)
            integer, intent(out)              :: nmerge

            if (.not. allocated(nppc)) then
                allocate(nppc(ncell))
                allocate(kc1(ncell))
                allocate(kc2(ncell))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc=0 !nppc(ic) will contain the number of parcels in grid cell ic

            vmin = vcell / dble(parcel%vmin_fraction)
            nmerge = 0

            ! Bin parcels in cells:
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

            ! if no small parcels are found, exit the subroutine
            if (nmerge == 0) then
                return
            endif

            ! Find arrays kc1(ic) & kc2(ic) which indicate the parcels in grid cell ic
            ! through i = node(k), for k = kc1(ic),kc2(ic):
            kc1(1)=1
            do ic=1,ncell-1
                kc1(ic+1)=kc1(ic)+nppc(ic)
            enddo

            !---------------------------------------------------------------------

            ! allocate arrays
            allocate(isma(0:nmerge + 1))
            allocate(ibig(nmerge + 1))

            isma = 0
            ibig = 0

            nmerge = 0

            kc2=kc1-1
            do i=1,n_parcels
                ic=loca(i)
                k=kc2(ic)+1
                node(k)=i
                kc2(ic)=k

                ! Form list of small parcels:
                if (parcels%volume(i) < vmin) then
                    nmerge = nmerge + 1
                    isma(nmerge) = i
                endif
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:
            j = 0
            ! j counts the actual number of mergers found
            do m=1, nmerge
                i0 = isma(m)
                ! Make sure parcel i0 is not already merging with a smaller one:
                x_small=parcels%position(i0,1)
                z_small=parcels%position(i0,2)
                ! Parcel i0 is small and should be merged; find closest other:
                ix0=mod(nint(dxi(1)*(x_small-lower(1))),nx) ! ranges from 0 to nx-1
                iz0=nint(dxi(2)*(z_small-lower(2)))         ! ranges from 0 to nz

                ! Grid point (ix0,iz0) is closest to parcel i0

                ! Initialise scaled squared distance between parcels and parcel index:
                ! In the loop below we want to minimise dsq/vmerge
                ! By storing dsqmin and vmergemin separately, we can avoid a division
                ! In the calculation we also want to ensure
                !   dsq/(a*b) < lambda_max/2
                ! This will ensure a merged parcel does not split again
                ! Since vmerge=pi*a*b, this implies
                !   dsq*pi < 0.5*parcel%lambda_max*vmerge
                ! This is ensured by initialising the minimisation
                ! with the values below
                ! Might seem a bit radical to take a large vmergemin and small dsqmin
                ! but computationally it is easy
                dsqmin = product(extent) !f12*parcel%lambda_max
!                 vmergemin=pi
                imin=0

                ! Loop over 8 cells surrounding (ix0,iz0):
                do iz=max(0,iz0-1),min(nz-1,iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do ix=ix0-1,ix0
                        ! Cell index (accounting for x periodicity):
                        ic=1+mod(nx+ix,nx)+nx*iz
                        ! Search nearby parcels for closest bigger one:
                        do k=kc1(ic),kc2(ic)
                            i=node(k)
                            if (parcels%volume(i) >= parcels%volume(i0)) then
                                ! we need to exclude self-merging due to v(i) >= v(i0)
                                if (i .ne. i0) then
                                    delz=parcels%position(i,2)-z_small
!                                     ! Avoid merger with another small parcel
!                                     vmerge=parcels%volume(i)+parcels%volume(i0) ! Summed area fraction:
!                                     ! Minimise dsq/vmerge
!                                     ! Prevent division in comparisons here by storing both
!                                     ! vmergemin and dsqmin
!                                     if (delz*delz < dsqmin) then
                                        ! works across periodic edge
                                        delx = get_delx(parcels%position(i,1), x_small)
                                        dsq=delz*delz+delx*delx
                                        if (dsq < dsqmin) then
                                            dsqmin=dsq
!                                             vmergemin=vmerge
                                            imin=i
                                        endif
!                                     endif
                                endif
                            endif
                        enddo
                    enddo
                enddo
                ! Store the index of the parcel found (or store 0 if none found):
                ibig(m) = imin
            enddo

            ! Pack isma and ibig such that we only have valid mergers
            j = 0
            do m = 1, nmerge
                if (ibig(m) > 0) then
                    ! (isma(m),ibig(m)) = (is,ib) is a valid merging pair
                    j = j + 1
                    is = isma(m)
                    isma(j) = is
                    ib = ibig(m)
                    ibig(j) = ib
                    ! Reuse the node array to identify incoming links
                    ! Reset here
                    node(is) = 0
                    node(ib) = 0
                endif
            enddo

            nmerge = j

            ! Identify number, type, identity of incoming mergers
            ! node(ib)==0 -> no incoming merger
            ! node(ib)==-1 -> smaller incoming merger, or multiple incoming mergers
            ! node(ib)>0 -> 1 incoming merger of equal size, keep track of id
            ! No re-packing here
            do m = 1, nmerge
                is = isma(m)
                ib = ibig(m)
                ! Reuse the node array to identify incoming links
                ! Keep  the number of incoming links
                if (parcels%volume(ib) > parcels%volume(is)) then
                    ! Big parcel has smaller incoming parcel
                    ! Exclude from merging
                    node(ib) = -1
                else if (node(ib) > 0) then
                    ! Big parcel has multiple equal size incoming parcels
                    ! Exclude from merging as at least 1 won't point back
                    node(ib) = -1
                else if (node(ib) == 0) then
                    ! Use node to identify identity of incoming parcel
                    ! IF THERE IS ONLY 1 SUCH PARCEL
                    node(ib) = is
                endif
            enddo

            ! Re-pack isma and ibig such that we only have valid mergers
            j = 0
            do m = 1, nmerge
                is = isma(m)
                if (node(is) == 0) then
                    ! No incoming mergers for is, keep merge
                    j = j + 1
                    isma(j) = is
                    ibig(j) = ibig(m)
                else if (node(is) > 0) then
                    ! Single incoming equal size mergers for this parcel
                    ib = ibig(m)
                    if (node(is) == ib) then
                        ! identity of incoming parcel is the "big" parcel
                        ! go ahead and merge
                        j = j + 1
                        isma(j) = is
                        ibig(j) = ib
                        ! but make sure the other parcel does not try to merge back
                        ! in case no other parcels point at it
                        node(ib) = -1
                    endif
                endif
                ! OTHERWISE, DO NOTHING
            enddo

            nmerge = j

!             do m = 1, nmerge
!                 print *, m, isma(m), ibig(m), node(isma(m)), node(ibig(m))
!             enddo

!             stop

        end subroutine find_nearest

end module parcel_nearest
