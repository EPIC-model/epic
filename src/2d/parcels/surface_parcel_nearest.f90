!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module surface_parcel_nearest
    use parcel_ops, only : get_delx
    use surface_parcel_container, only : surface_parcel_container_type
    use parameters, only : dx, dxi, lower, extent, nx, lmin, max_num_surf_parcels
    implicit none

    integer:: merge_nearest_timer, merge_tree_resolve_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:),kc2(:)
    integer, allocatable :: loca(:)
    integer, allocatable :: node(:)

    ! Logicals used to determine which mergers are executed
    ! Integers above could be reused for this, but this would
    ! make the algorithm less readable
    logical, allocatable :: l_leaf(:)
    logical, allocatable :: l_available(:)
    logical, allocatable :: l_first_merged(:) ! indicates parcels merged in first stage

#ifndef NDEBUG
    ! Logicals that are only needed for sanity checks
    logical, allocatable :: l_merged(:)! SANITY CHECK ONLY
    logical, allocatable :: l_small(:) ! SANITY CHECK ONLY
    logical, allocatable :: l_close(:) ! SANITY CHECK ONLY
#endif

    logical :: l_continue_iteration, l_do_merge

    !Other variables:
    double precision:: delx, dsq, dsqmin, x_small
    integer:: ic, is, ij, k, m, j, n
    integer:: ix, ix0

    public :: find_nearest, merge_nearest_timer, merge_tree_resolve_timer

    contains

        ! @param[out] isma indices of small parcels
        ! @param[out] iclo indices of close parcels
        ! @param[out] nmerge the array size of isma and iclo
        ! @post
        !   - isma must be sorted in ascending order
        !   - isma and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(n_par, spar, isma, iclo, nmerge)
            integer,                             intent(in)  :: n_par
            type(surface_parcel_container_type), intent(in)  :: spar
            integer, allocatable,                intent(out) :: isma(:)
            integer, allocatable,                intent(out) :: iclo(:)
            integer,                             intent(out) :: nmerge

            if (.not. allocated(nppc)) then
                allocate(nppc(nx))
                allocate(kc1(nx))
                allocate(kc2(nx))
                allocate(loca(max_num_surf_parcels))
                allocate(node(max_num_surf_parcels))
                allocate(l_leaf(max_num_surf_parcels))
                allocate(l_available(max_num_surf_parcels))
                allocate(l_first_merged(max_num_surf_parcels))
#ifndef NDEBUG
                allocate(l_merged(max_num_surf_parcels))
                allocate(l_small(max_num_surf_parcels))
                allocate(l_close(max_num_surf_parcels))
#endif
            endif

            nmerge = 0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(ij) will contain the number of parcels in grid cell ij

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_par
                ix = mod(int(dxi(1) * (spar%position(n) - lower(1))), nx)

                ! Cell index of parcel:
                ij = 1 + ix !This runs from 1 to nx

                ! Accumulate number of parcels in this grid cell:
                nppc(ij) = nppc(ij) + 1

                ! Store grid cell that this parcel is in:
                loca(n) = ij

                if (spar%length(n) < lmin) then
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
            do ij = 1, nx-1
                kc1(ij+1) = kc1(ij) + nppc(ij)
            enddo

            kc2 = kc1 - 1
            j = 0
            do n = 1, n_par
                ij = loca(n)
                k = kc2(ij) + 1
                node(k) = n
                kc2(ij) = k

                if (spar%length(n) < lmin) then
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

            ! SB: do not use temporary (j) index here, so we will be able to parallelise.
            ! Rather, stop if no nearest parcel found  in surrounding grid boxes
            do m = 1, nmerge
                is = isma(m)
                x_small = spar%position(is)
                ! Parcel "is" is small and should be merged; find closest other:
                ix0 = mod(nint(dxi(1) * (x_small - lower(1))), nx) ! ranges from 0 to nx-1

                ! Grid point (ix0,iz0) is closest to parcel "is"

                dsqmin = extent(1)
                ic = 0

                ! Loop over 2 cells surrounding (ix0):
                do ix = ix0-1, ix0
                    ! Cell index (accounting for x periodicity):
                    ij = 1 + mod(nx + ix, nx)
                    ! Search small parcels for closest other:
                    do k = kc1(ij), kc2(ij)
                        n = node(k)
                        if (n .ne. is) then
                            delx = get_delx(spar%position(n), x_small) ! works across periodic edge
                            ! Minimise dsqmin
                            dsq = delx * delx
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
                end if

                ! Store the index of the parcel to be potentially merged with:
                isma(m) = is
                iclo(m) = ic
                l_first_merged(is)=.false.
                l_first_merged(ic)=.false.
            enddo

#ifndef NDEBUG
            write(*,*) 'start merging, nmerge='
            write(*,*) nmerge
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
                ! only consider links that still may be merging
                ! reset relevant properties
                if(.not. l_first_merged(is)) then
                  ic = iclo(m)
                  l_leaf(is)=.true.
                  l_available(ic)=.true.
                end if
              end do
              ! determine leaf parcels
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_first_merged(is)) then
                  ic = iclo(m)
                  l_leaf(ic)=.false.
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! filter out parcels that are "unavailable" for merging
              do m = 1, nmerge
                is = isma(m)
                if(.not. l_first_merged(is)) then
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
                if(.not. l_first_merged(is)) then
                  ic = iclo(m)
                  if(l_leaf(is) .and. l_available(ic)) then
                    l_continue_iteration=.true. ! merger means continue iteration
                    l_first_merged(is)=.true.
                    l_first_merged(ic)=.true.
                  end if
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
            end do

            ! Second stage, related to dual links
            do m = 1, nmerge
              is = isma(m)
              if(.not. l_first_merged(is)) then
                if(l_leaf(is)) then ! set in last iteration of stage 1
                  ic = iclo(m)
                  l_available(ic)=.true.
                end if
              end if
            end do

            ! PARALLEL: COMMUNICATION WILL BE NEEDED HERE

            ! Second stage (hard to parallelise with openmp?)
            j=0
            do m = 1, nmerge
              is = isma(m)
              ic = iclo(m)
              l_do_merge=.false.
              if(l_first_merged(is) .and. l_leaf(is)) then
                ! previously identified mergers: keep
                l_do_merge=.true.
#ifndef NDEBUG
                ! sanity check on first stage mergers
                ! parcel cannot be both initiator and receiver in stage 1
                if(l_leaf(ic)) then
                  write(*,*) 'first stage error'
                endif
#endif
              elseif(.not. l_first_merged(is)) then
                if(l_leaf(is)) then
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
                endif
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
            write(*,*) 'after second stage, nmerge='
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
            do n = 1, n_par
              if (spar%length(n) < lmin) then
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

end module surface_parcel_nearest
