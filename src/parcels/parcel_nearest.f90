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
                            if (parcels%volume(i) >= parcels%volume(i0)) then
                                ! we need to exclude self-merging due to v(i) >= v(i0)
                                if (i .ne. i0) then
                                    delz=parcels%position(i,2)-z_small
                                    ! Avoid merger with another small parcel
                                    vmerge=parcels%volume(i)+parcels%volume(i0) ! Summed area fraction:
                                    ! Minimise dsq/vmerge
                                    ! Prevent division in comparisons here by storing both
                                    ! vmergemin and dsqmin
                                    if (delz*delz*vmergemin < dsqmin*vmerge) then
                                        ! works across periodic edge
                                        delx = get_delx(parcels%position(i,1), x_small)
                                        dsq=delz*delz+delx*delx
                                        if (dsq*vmergemin < dsqmin*vmerge) then
                                            dsqmin=dsq
                                            vmergemin=vmerge
                                            imin=i
                                        endif
                                    endif
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
                    ibig(j) = ibig(m)
                    ! Reuse the node array to identify incoming links
                    ! Reset here
                    node(is) = 0
                endif
            enddo

            !---------------------------------------------------------------------
            ! Link merging parcels (isma(m),ibig(m)) = (is,ib) by node(is) = ib:
            do m = 1, nmerge
                node(isma(m)) = ibig(m) ! if ibig(m) = 0, parcel isma(m) is isolated
            enddo
            ! This re-uses the node array above - no extra memory is required.

            ! First, eliminate mergers between equal-sized parcels when they
            ! do NOT point to each other:
            do m = 1, nmerge
                is = isma(m)
                ib = node(is)
                if (ib > 0) then
                    if (v(is) == v(ib)) then
                        if (node(ib) .ne. is) then
                            node(is) = 0
                        endif
                    endif
                endif
            enddo
            ! An example here is A <-> B <- C, i.e. node(A) = B, node(B) = A and
            ! node(C) = B.  The above eliminates C, retaining A <-> B.

            ! Next, forbid a strictly larger parcel "ib" from ever merging with
            ! another of equal or larger size; do this by setting node(ib) = 0:
            do m = 1, nmerge
                ib = ibig(m)
                if (ib > 0) then
                    ! If larger parcel in a pair; forbid from further mergers:
                    if (v(ib) > v(isma(m))) then
                        node(ib) = 0
                    endif
                endif
            enddo
            ! The above two loops do not depend on the order parcels are found
            ! in the list over m.  Any order eliminates bigger parcels in a chain
            ! like A -> B -> C -> D ... (C & D etc are eliminated, leaving A -> B).

            ! Pack isma and ibig such that we only have valid mergers:
            j = 0
            do m = 1, nmerge
                is = isma(m)
                ib = node(is)
                if (ib > 0) then
                    j = j + 1
                    isma(j) = is
                    ibig(j) = ib
                    if (node(ib) == is) then
                        node(ib) = 0
                    endif
                    ! This prevents the merger between equal-sized parcels from being
                    ! counted twice; ONLY when v(is) = v(ib) would node(ib) = is;
                    ! resetting node(ib) = 0 ensures this case is only counted once.
                    ! That is, in A <-> B, one of the links is broken, leaving either
                    ! A -> B or B -> A; but this symmetric as other parcel is involved.
                endif
            enddo

            ! Actual total number of valid binary mergers:
            nmerge = j

!             do m = 1, nmerge
!                 print *, m, isma(m), ibig(m), node(isma(m)), node(ibig(m))
!             enddo

!             stop

        end subroutine find_nearest

end module parcel_nearest
