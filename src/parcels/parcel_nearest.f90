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
    logical :: l_leaf(max_num_parcels)
    logical :: l_available(max_num_parcels)
    logical :: l_first_stage(max_num_parcels)

#ifndef NDEBUG
    logical :: l_merged(max_num_parcels)! SANITY CHECK ONLY
    logical :: l_small(max_num_parcels) ! SANITY CHECK ONLY
    logical :: l_close(max_num_parcels) ! SANITY CHECK ONLY
#endif

    logical :: l_continue_iteration, l_do_merge

    ! Logicals that are only needed for sanity checks

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
#ifndef NDEBUG
                l_merged(n)=.false.! SANITY CHECK ONLY
                l_small(n)=.false. ! SANITY CHECK ONLY
                l_close(n)=.false. ! SANITY CHECK ONLY
#endif
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
                  l_first_stage(ic)=.false.
                  l_first_stage(is)=.false.
                endif

            enddo


            ! Actual total number of mergers:
            nmerge = j
#ifndef NDEBUG
            write(*,*) 'start merging, nmerge='
            write(*,*) nmerge
            write(*,*) 'first stage, nmerge='
#endif

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
                if(.not. l_first_stage(is)) then
                  ic = iclo(m)
                  l_leaf(is)=.true.
                  l_available(ic)=.true.
                end if
              end do
              ! determine leaf parcels
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_first_stage(is)) then
                  ic = iclo(m)
                  l_leaf(ic)=.false.
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! filter out for big parcels with height>1
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_first_stage(is)) then
                   if(.not. l_leaf(is)) then
                     ic = iclo(m)
                     l_available(ic)=.false.
                   end if
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! identify mergers in this iteration
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_first_stage(is)) then
                  ic = iclo(m)
                  if(l_leaf(is) .and. l_available(ic)) then
                    l_continue_iteration=.true. ! merger means continue iteration
                    l_first_stage(ic)=.true.
                    l_first_stage(is)=.true.
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
                if(.not.(l_first_stage(is) .and. (.not. l_first_stage(ic)))) then
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic
                endif
              enddo
              nmerge = j
#ifndef NDEBUG
              write(*,*) nmerge
#endif
            end do

            ! Second stage, related to dual links
            do m = 1, nmerge
              is = isma(m)
              if(.not. l_first_stage(is)) then
                if(l_leaf(is)) then
                  ic = iclo(m)
                  l_available(ic)=.true.
                end if
              end if
            end do

            ! PARALLEL: COMMUNICATION WILL BE NEEDED HERE

            j=0
            do m = 1, nmerge
              is = isma(m)
              ic = iclo(m)
              l_do_merge=.false.
              if(l_first_stage(is) .and. l_first_stage(ic)) then
                 ! previously identified mergers: keep
                 l_do_merge=.true.
              elseif(l_leaf(is)) then
                 ! links from leafs
                 l_do_merge=.true.
              elseif(.not. l_available(is)) then
                 ! Above means parcels that have been made 'available' do not keep outgoing links
                 if(l_available(ic)) then
                   ! merge this parcel into ic along with the leaf parcels
                   l_do_merge=.true.
                 else
                   ! isolated dual link
                   ! Don't keep current link
                   ! But make small parcel available so other parcel can merge with it
                   ! THIS NEEDS THINKING ABOUT A PARALLEL IMPLEMENTATION
                   ! This could be based on the other parcel being outside the domain
                   ! And a "processor order"
                   l_available(is)=.true.
                 endif
#ifndef NDEBUG
              elseif(l_first_stage(is)) then
                 write(*,*) 'first stage error'
              elseif(l_first_stage(ic)) then
                 write(*,*) 'second stage error'
#endif
              endif

              if(l_do_merge) then
                   j = j + 1
                   isma(j) = is
                   iclo(j) = ic
#ifndef NDEBUG
                   l_merged(is)=.true.
                   l_merged(ic)=.true.
                   l_small(is)=.true.
                   l_close(ic)=.true.
#endif
              end if
            end do
            nmerge = j

#ifndef NDEBUG
            write(*,*) 'second stage'
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
              if(.not. l_merged(isma(m))) write(*,*) 'merge_error: isma(m) not merged, m=', m
              if(.not. l_merged(iclo(m))) write(*,*) 'merge_error: iclo(m) not merged, m=', m
              if(.not. l_small(isma(m))) write(*,*) 'merge_error: isma(m) not marked as small, m=', m
              if(.not. l_close(iclo(m))) write(*,*) 'merge_error: iclo(m) not marked as close, m=', m
              if(l_close(isma(m))) write(*,*) 'merge_error: isma(m) both small and close, m=', m
              if(l_small(iclo(m))) write(*,*) 'merge_error: iclo(m) both small and close, m=', m
            end do

            ! 2. CHECK MERGING PARCELS
            do n = 1, n_parcels
              if (parcels%volume(n) < vmin) then
                if(.not. l_merged(n)) write(*,*) 'merge_error: parcel n not merged (should be), n=', n
                if(.not. (l_small(n) .or. l_close(n))) write(*,*) 'merge_error: parcel n not small or close (should be), n=', n
                if(l_small(n) .and. l_close(n)) write(*,*) 'merge_error: parcel n both small and close, n=', n
              else
                if(l_small(n)) write(*,*) 'merge_error: parcel n small (should not be), n=', n
                if(l_merged(n) .and. (.not. l_close(n))) write(*,*) 'merge_error: parcel n merged (should not be), n=', n
              end if
            enddo
#endif

        end subroutine find_nearest

end module parcel_nearest
