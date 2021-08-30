!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : pi, f12, max_num_parcels
    use parcel_container, only : parcels, n_parcels, get_delx
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, nz, vmin
    use options, only : parcel
    use merge_sort

    implicit none

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:),kc2(:)
    integer :: loca(max_num_parcels)
    integer :: node(max_num_parcels)

    ! Logicals used to determine which mergers are executed
    ! Integers above could be reused for this, but this would
    ! make the algorithm less readable
    logical :: l_is_leaf(max_num_parcels)
    logical :: l_is_available(max_num_parcels)
    logical :: l_is_merged(max_num_parcels)
    logical :: l_first_stage(max_num_parcels)
    logical :: l_continue_iteration

    ! Logicals that are only needed for sanity checks
    logical :: l_is_small(max_num_parcels) ! SANITY CHECK ONLY
    logical :: l_is_close(max_num_parcels) ! SANITY CHECK ONLY

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

            ! Find arrays kc1(ij) & kc2(ij) which indicate the parcels in grid cell ij
            ! through n = node(k), for k = kc1(ij),kc2(ij):
            kc1(1) = 1
            do ij = 1, ncell-1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            j = 0
            do n = 1, n_parcels
                ij = loca(n)
                k = kc2(ij) + 1
                node(k) = n
                kc2(ij) = k

                if (parcels%volume(n) < vmin) then
                    j = j + 1
                    isma(j) = n
                endif
                l_is_close(n)=.false. ! SANITY CHECK ONLY
                l_is_small(n)=.false. ! SANITY CHECK ONLY
                l_is_merged(n)=.false.! SANITY CHECK ONLY (OTHERWISE COULD BE SET IN LOOP BELOW)
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:
            j=0
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
                if (ic>0) then
                  j=j+1
                  isma(j) = is
                  iclo(j) = ic
                  l_is_merged(is)=.false.
                  l_is_merged(ic)=.false.
                  l_first_stage(ic)=.false.
                  l_first_stage(is)=.false.
                endif

            enddo

            write(*,*) 'start'

            ! Actual total number of mergers:
            nmerge = j
            write(*,*) nmerge

            ! First implementation of iterative algorithm

            ! Could be improved by keeping track of "finalised mergers"
            ! But going for more readable implementation here first.

            ! First, iterative, stage
            l_continue_iteration=.true.
            do while(l_continue_iteration)
              l_continue_iteration=.false.
              ! reset relevant properties for candidate mergers
              do m = 1, nmerge
                is = isma(m)
                ! only consider links that need merging
                ! reset relevant properties
                if(.not. l_is_merged(is)) then
                  ic = iclo(m)
                  l_is_leaf(is)=.true.
                  l_is_available(ic)=.true.
                end if
              end do
              ! determine leaf parcels
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_is_merged(is)) then
                  ic = iclo(m)
                  l_is_leaf(ic)=.false.
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! filter out for big parcels with height>1
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_is_merged(is)) then
                   if(.not. l_is_leaf(is)) then
                     ic = iclo(m)
                     l_is_available(ic)=.false.
                   end if
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! identify mergers in this iteration
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_is_merged(is)) then
                  ic = iclo(m)
                  if(l_is_leaf(is) .and. l_is_available(ic)) then
                    l_continue_iteration=.true. ! merger means continue iteration
                    l_is_merged(is)=.true.
                    l_is_merged(ic)=.true. ! note multiple parcels can merge into ic
                    l_first_stage(ic)=.true.
                    l_first_stage(is)=.true.
                    l_is_small(is)=.true. !SANITY CHECK ONLY
                    l_is_close(ic)=.true. !SANITY CHECK ONLY
                  end if
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! remove eliminated links
              j = 0
              do m = 1, nmerge
                is = isma(m)
                ic = iclo(m)
                ! keep the link, unless small parcel merges but close one does not
                if(.not.(l_is_merged(is) .and. (.not. l_is_merged(ic)))) then
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic
                endif
              enddo
              nmerge = j
              write(*,*) nmerge
            end do

            write(*,*) 'second stage'

            ! Second stage, related to dual links
            do m = 1, nmerge
              is = isma(m)
              if(.not. l_is_merged(is)) then
                if(l_is_leaf(is)) then
                  ic = iclo(m)
                  l_is_available(ic)=.true.
                end if
              end if
            end do

            ! PARALLEL: COMMUNICATION WILL BE NEEDED HERE

            j=0
            do m = 1, nmerge
              is = isma(m)
              ic = iclo(m)
              if(l_first_stage(is) .and. l_first_stage(ic)) then
                 ! previously identified mergers: keep
                 j = j + 1
                 isma(j) = is
                 iclo(j) = ic
              elseif(l_is_leaf(is)) then
                 ! links from leafs
                 j = j + 1
                 isma(j) = is
                 iclo(j) = ic
                 l_is_merged(is)=.true.
                 l_is_merged(ic)=.true.
                 l_is_small(is)=.true. !SANITY CHECK ONLY
                 l_is_close(ic)=.true. !SANITY CHECK ONLY
              elseif(.not. l_is_available(is)) then
                 if(l_is_available(ic)) then
                   ! merge this parcel into ic along with the leaf parcels
                   j = j + 1
                   isma(j) = is
                   iclo(j) = ic
                   l_is_merged(is)=.true.
                   l_is_merged(ic)=.true.
                   l_is_small(is)=.true. !SANITY CHECK ONLY
                   l_is_close(ic)=.true. !SANITY CHECK ONLY
                 else
                   ! isolated dual link
                   ! Don't keep current link
                   ! But make small parcel available so other parcel can merge with it
                   ! THIS NEEDS THINKING ABOUT A PARALLEL IMPLEMENTATION
                   ! This could be based on the other parcel being outside the domain
                   ! And a "processor order"
                   l_is_available(is)=.true.
                 endif
              elseif(l_first_stage(is)) then
                 write(*,*) 'first stage error'
              elseif(l_first_stage(ic)) then
                 write(*,*) 'second stage error'
              endif
              ! This means parcels that have been made 'available' do not keep outgoing links
            end do
            nmerge = j

            write(*,*) nmerge

            write(*,*) 'finished'

            ! MORE SANITY CHECKS
            ! CHECK ISMA ORDER
            do m = 1, nmerge
              if(.not. (isma(m)>isma(m-1))) then
                write(*,*) 'isma order broken'
              end if
            end do

            ! 1. CHECK RESULTING MERGERS
            do m = 1, nmerge
              if(.not. l_is_merged(isma(m))) write(*,*) 'merge_error: isma(m) not merged', m
              if(.not. l_is_merged(iclo(m))) write(*,*) 'merge_error: iclo(m) not merged', m
              if(.not. l_is_small(isma(m))) write(*,*) 'merge_error: isma(m) not marked as small', m
              if(.not. l_is_close(iclo(m))) write(*,*) 'merge_error: iclo(m) not marked as close', m
              if(l_is_close(isma(m))) write(*,*) 'merge_error: isma(m) both small and close', m
              if(l_is_small(iclo(m))) write(*,*) 'merge_error: iclo(m) both small and close', m
            end do

            ! 2. CHECK MERGING PARCELS
            do n = 1, n_parcels
              if (parcels%volume(n) < vmin) then
                if(.not. l_is_merged(n)) write(*,*) 'merge_error: parcel n not merged (should be)', n
                if(.not. (l_is_small(n) .or. l_is_close(n))) write(*,*) 'merge_error: parcel n not small or close (should be)', n
                if(l_is_small(n) .and. l_is_close(n)) write(*,*) 'merge_error: parcel n both small and close', n
              else
                if(l_is_small(n)) write(*,*) 'merge_error: parcel n small (should not be)', n
                if(l_is_merged(n) .and. (.not. l_is_close(n))) write(*,*) 'merge_error: parcel n merged (should not be)', n
              end if
            enddo


        end subroutine find_nearest

end module parcel_nearest
