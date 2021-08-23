!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : pi, f12, max_num_parcels
    use parcel_container, only : parcels, n_parcels, get_delx
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, nz, vmin
    use options, only : parcel

    implicit none

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:),kc2(:)
    integer :: loca(max_num_parcels)
    integer :: node(max_num_parcels)

    !Other variables:
    double precision:: delx,delz,dsq,dsqmin,x_small,z_small
    integer:: ic,is,ij,k,m,j, n
    integer:: ix,iz,ix0,iz0

    public :: find_nearest

    contains

        subroutine find_nearest(isma, iclo, nmerge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer, intent(out) :: nmerge

            if (.not. allocated(nppc)) then
                allocate(nppc(ncell))
                allocate(kc1(ncell))
                allocate(kc2(ncell))
            endif

            nmerge = 0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(ij) will contain the number of parcels in grid cell ij

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels
                ix = int(dxi(1) * (parcels%position(n,1) - lower(1)))
                iz = int(dxi(2) * (parcels%position(n,2) - lower(2)))

                ! Cell index of parcel:
                ij = 1 + ix + nx * iz !This runs from 1 to ncell

                ! Accumulate number of parcels in this grid cell:
                nppc(ij) = nppc(ij) + 1

                ! Store grid cell that this parcel is in:
                loca(n) = ij

                if (parcels%volume(n) < vmin) then
                    nmerge = nmerge + 1
                endif
            enddo

            if (nmerge == 0) then
                return
            endif

            ! allocate arrays
            allocate(isma(0:nmerge))
            allocate(iclo(nmerge))

            isma = 0
            iclo = 0
            nmerge = 0

            ! Find arrays kc1(ij) & kc2(ij) which indicate the parcels in grid cell ij
            ! through n = node(k), for k = kc1(ij),kc2(ij):
            kc1(1) = 1
            do ij = 1, ncell-1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            do n = 1, n_parcels
                ij = loca(n)
                k = kc2(ij) + 1
                node(k) = n
                kc2(ij) = k

                if (parcels%volume(n) < vmin) then
                    nmerge = nmerge + 1
                    isma(nmerge) = n
                endif
                ! reset loca (for use later)
                loca(n) = -1
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:
            j = 0
            ! j counts the actual number of mergers found
            do m = 1, nmerge
                is = isma(m)
                x_small = parcels%position(is, 1)
                z_small = parcels%position(is, 2)
                ! Parcel "is" is small and should be merged; find closest other:
                ix0 = mod(nint(dxi(1) * (x_small - lower(1))), nx) ! ranges from 0 to nx-1
                iz0 = nint(dxi(2) * (z_small - lower(2)))          ! ranges from 0 to nz

                ! Grid point (ix0,iz0) is closest to parcel "is"

                dsqmin = product(extent)
                ic = 0

                ! Loop over 8 cells surrounding (ix0,iz0):
                do iz = max(0, iz0-1), min(nz-1, iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do ix = ix0-1, ix0
                        ! Cell index (accounting for x periodicity):
                        ij = 1 + mod(nx + ix, nx) + nx * iz
                        ! Search small parcels for closest other:
                        do k = kc1(ij), kc2(ij)
                            n = node(k)
                            if (n .ne. is) then
                                delz = parcels%position(n,2) - z_small
                                delx = get_delx(parcels%position(n,1), x_small) ! works across periodic edge
                                ! Minimise dsqmin
                                dsq = delz * delz + delx * delx
                                if (dsq < dsqmin) then
                                    dsqmin = dsq
                                    ic = n
                                endif
                            endif
                        enddo
                    enddo
                enddo

                ! Store the index of the parcel to be potentially merged with:
                j = j + 1
                isma(j) = is
                iclo(j) = ic
                loca(is) = ic
            enddo
            ! Actual total number of mergers:
            nmerge = j

            j = 0
            do m = 1, nmerge
                ic = iclo(m)
                is = isma(m)

                if (loca(ic) == -1) then
                    ! ic is only big
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic
                else if (loca(ic) == is) then
                    ! dual link or double bond
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic

                    loca(is) = -2
                    loca(ic) = -2
                ! else
                    ! "ic" is small but not a dual link; "is" is not merged
                endif

            enddo
            nmerge = j

        end subroutine find_nearest

end module parcel_nearest
