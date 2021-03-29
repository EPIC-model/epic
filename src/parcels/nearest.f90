!===================================================================
!         Finds the parcels nearest every "small" parcel

!         Initialise data first using ranpar.f90
!===================================================================

module nearest
    use constants, only : pi, max_num_parcels
    use parcel_container, only : parcels, n_parcels
    use constants, only : pi
    use parameters
    use fields, only : get_mesh_spacing

    implicit none

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:),kc2(:) !ncell = nx*nz
    integer :: loc(max_num_parcels)
    integer :: node(max_num_parcels)
!     integer :: isma(max_num_parcels/8)
!     integer :: ibig(max_num_parcels/8)
!     integer :: nmerge
    logical :: merge(max_num_parcels)

    !Other variables:
    double precision:: vmin, delx,delz,dsq,dscmax,dscmin,vmerge,vmergemin,x_store,z_store
    integer:: i,ic,i0,imin,k,m
    integer:: ix,iz,ix0,iz0

    contains

        ! if ibig(n) is zero, no parcel found for isma(n)
        subroutine find_nearest(isma, ibig, nmerge)
            integer, intent(out) :: isma(max_num_parcels)
            integer, intent(out) :: ibig(max_num_parcels)
            integer, intent(out) :: nmerge
            integer          :: nx, nz, ncell
            double precision :: dx(2), dxi(2)
            !Inverse of domain half width:
            double precision :: hlxi

            hlxi = 0.5 * mesh%extent(1)

            nx = mesh%grid(1)
            nz = mesh%grid(2)
            ncell = nx * nz

            if (.not. allocated(nppc)) then
                allocate(nppc(ncell))
                allocate(kc1(ncell))
                allocate(kc2(ncell))
            endif

            dx = get_mesh_spacing()
            dxi = 1.0 / dx

            ! maximum squared distance
            dscmax = 0.5 * parcel_info%lambda

            vmin = product(dx) / (parcel_info%n_per_cell * parcel_info%max_splits)

            ! These parcels are marked for merger:
            merge(1:n_parcels)=(parcels%volume(1:n_parcels, 1) < vmin)
            nmerge=0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc=0 !nppc(ic) will contain the number of parcels in grid cell ic

            ! Form list of small parcels:
            do i=1,n_parcels
                if (merge(i)) then
                    nmerge=nmerge+1
                    isma(nmerge)=i
                else
                    ix=int(dxi(1)*(parcels%position(i,1)-mesh%origin(1)))
                    iz=int(dxi(2)*(parcels%position(i,2)-mesh%origin(2)))

                    ! Cell index of parcel:
                    ic=1+ix+nx*iz !This runs from 1 to ncell

                    ! Accumulate number of parcels in this grid cell:
                    nppc(ic)=nppc(ic)+1

                    ! Store grid cell that this parcel is in:
                    loc(i)=ic
                endif
            enddo

            ! Find arrays kc1(ic) & kc2(ic) which indicate the parcels in grid cell ic
            ! through i = node(k), for k = kc1(ic),kc2(ic):
            kc1(1)=1
            do ic=1,ncell-1
                kc1(ic+1)=kc1(ic)+nppc(ic)
            enddo

            kc2=kc1-1
            do i=1,n_parcels
                if (.not. merge(i)) then
                    ic=loc(i)
                    k=kc2(ic)+1
                    node(k)=i
                    kc2(ic)=k
                end if
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:
            do m=1,nmerge
                i0=isma(m)
                x_store=parcels%position(i0,1)
                z_store=parcels%position(i0,2)
                ! Parcel i0 is small and should be merged; find closest other:
                ix0=mod(nint(dxi(1)*(x_store-mesh%origin(1))),nx) ! ranges from 0 to nx-1
                iz0=nint(dxi(2)*(z_store-mesh%origin(2)))         ! ranges from 0 to nz
                ! Grid point (ix0,iz0) is closest to parcel i0

                ! Initialise scaled squared distance between parcels and parcel index:
                ! ensure dsq*pi < 0.5*parcel_info%lambda*vmerge
                ! Might seem a bit radical to take a large vmergemin and small dscmin
                ! but computationally easy
                dscmin=0.5*parcel_info%lambda
                vmergemin=pi
                imin=0

                ! Loop over 8 cells surrounding (ix0,iz0):
                do iz=max(0,iz0-1),min(nz-1,iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do ix=ix0-1,ix0
                        ! Cell index (accounting for x periodicity):
                        ic=1+mod(nx+ix,nx)+nx*iz
                        ! Search parcels for closest:
                        do k=kc1(ic),kc2(ic)
                            i=node(k)
                            delz=parcels%position(i,2)-z_store
                            ! Avoid merger with another small parcel
                            vmerge=parcels%volume(i, 1)+parcels%volume(i0, 1) ! Summed area fraction:
                            ! Prevent division in all comparisons here
                            if (delz*delz*vmergemin < dscmin*vmerge) then
                                delx=parcels%position(i,1)-x_store
                                delx=delx-mesh%extent(1)*dble(int(delx*hlxi)) ! works across periodic edge
                                dsq=delz*delz+delx*delx
                                if (dsq*vmergemin < dscmin*vmerge) then
                                    dscmin=dsq
                                    vmergemin=vmerge
                                    imin=i
                                endif
                            endif
                        enddo
                    enddo
                enddo
                ! Store the index of the parcel to be merged with:
                ibig(m)=imin
            enddo
        end subroutine find_nearest

end module nearest
