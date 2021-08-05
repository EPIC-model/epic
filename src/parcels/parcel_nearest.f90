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
    double precision:: vmin, delx,delz,dsq,dsqmin,vmerge,vmergemin,x_small,z_small
    integer :: i,ic,i0,imin,k,m,j
    integer :: ix,iz,ix0,iz0

    public :: find_nearest

    contains

        ! if ibig(n) is zero, no parcel found for isma(n)
        subroutine find_nearest(isma, ibig, nmerge)
            integer, intent(out) :: isma(0:n_parcels)
            integer, intent(out) :: ibig(n_parcels)
            integer, intent(out) :: nmerge
            logical              :: avail(n_parcels) ! indicates that parcel i is available for merger

            if (.not. allocated(nppc)) then
                allocate(nppc(ncell))
                allocate(kc1(ncell))
                allocate(kc2(ncell))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc=0 !nppc(ic) will contain the number of parcels in grid cell ic

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
            enddo

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
            enddo

            !---------------------------------------------------------------------
            ! Form list of small parcels:
            nmerge = 0
            vmin = vcell / dble(parcel%vmin_fraction)

            do i = 1, n_parcels
                if (parcels%volume(i) < vmin) then
                    nmerge = nmerge + 1
                    isma(nmerge) = i
                    avail(i) = .true.
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
                if (avail(i0)) then
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
                    dsqmin=f12*parcel%lambda_max
                    vmergemin=pi
                    imin=0

                    ! Loop over 8 cells surrounding (ix0,iz0):
                    do iz=max(0,iz0-1),min(nz-1,iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                        do ix=ix0-1,ix0
                            ! Cell index (accounting for x periodicity):
                            ic=1+mod(nx+ix,nx)+nx*iz
                            ! Search nearby parcels for closest bigger one:
                            do k=kc1(ic),kc2(ic)
                                i=node(k)
                                if (parcels%volume(i) > parcels%volume(i0)) then
                                    delz=parcels%position(i,2)-z_small
                                    ! Avoid merger with another small parcel
                                    vmerge=parcels%volume(i)+parcels%volume(i0) ! Summed area fraction:
                                    ! Minimise dsq/vmerge
                                    ! Prevent division in comparisons here by storing both
                                    ! vmergemin and dsqmin
                                    if (delz*delz*vmergemin < dsqmin*vmerge) then
                                        delx = get_delx(parcels%position(i,1), x_small) ! works across periodic edge
                                        dsq=delz*delz+delx*delx
                                        if (dsq*vmergemin < dsqmin*vmerge) then
                                            dsqmin=dsq
                                            vmergemin=vmerge
                                            imin=i
                                        endif
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                    if (imin .gt. 0) then
                        ! Store the indices of the merging parcels:
                        j = j + 1
                        isma(j) = i0
                        ibig(j) = imin
                        ! Forbid the bigger parcel (imin) from merging:
                        avail(imin) = .false.
                    endif
                endif
            enddo
            ! Actual total number of binary mergers:
            nmerge = j

            ! Pack isma and ibig such that we only have valid mergers
            j = 0
            do m = 1, nmerge
                if (avail(isma(m))) then
                    j = j + 1
                    isma(j) = isma(m)
                    ibig(j) = ibig(m)
                endif
            enddo

            ! Total number of merges:
            nmerge = j
        end subroutine find_nearest

end module parcel_nearest
