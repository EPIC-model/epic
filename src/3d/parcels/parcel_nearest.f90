!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : pi, f12
    use parcel_container, only : parcels, n_parcels, get_delx, get_dely
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use mpi_layout, only : box
    use timer, only : start_timer, stop_timer

    implicit none

    integer:: merge_nearest_timer, merge_tree_resolve_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:), kc2(:)
    integer, allocatable :: lloca(:)
    integer, allocatable :: gloca(:)
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
    double precision:: delx, dely, delz, dsq, dsqmin, x_small, y_small, z_small
    integer :: ic, is, k, m, j, n
    integer :: lijk ! local ijk index
    integer :: gijk ! global ijk index
    integer :: ix, iy, iz, ix0, iy0, iz0
    integer :: n_halo_small(8)  ! number of small parcels in halo regions

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
        subroutine find_nearest(isma, iclo, nmerge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer, intent(out) :: nmerge

            call start_timer(merge_nearest_timer)

            if (.not. allocated(nppc)) then
                allocate(nppc(box%ncell))
                allocate(kc1(box%ncell))
                allocate(kc2(box%ncell))
                allocate(lloca(max_num_parcels))
                allocate(gloca(max_num_parcels))
                allocate(node(max_num_parcels))
                allocate(l_leaf(max_num_parcels))
                allocate(l_available(max_num_parcels))
                allocate(l_first_merged(max_num_parcels))
#ifndef NDEBUG
                allocate(l_merged(max_num_parcels))
                allocate(l_small(max_num_parcels))
                allocate(l_close(max_num_parcels))
#endif
            endif

            nmerge = 0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(lijk) will contain the number of parcels in grid cell gijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels
                ix = int(dxi(1) * (parcels%position(1, n) - lower(1)))
                iy = int(dxi(2) * (parcels%position(2, n) - lower(2)))
                iz = int(dxi(3) * (parcels%position(3, n) - lower(3)))

                ! Cell index of parcel:
                gijk = 1 + ix +     nx * iy +     nx *     ny * iz
                lijk = 1 + ix + box%nx * iy + box%nx * box%ny * iz !This runs from 1 to box%ncell

                ! Accumulate number of parcels in this grid cell:
                nppc(lijk) = nppc(lijk) + 1

                ! Store grid cell that this parcel is in:
                lloca(n) = lijk
                gloca(n) = gijk

                if (parcels%volume(n) < vmin) then
                    nmerge = nmerge + 1
                endif
            enddo

            call MPI_Allreduce(MPI_IN_PLACE, nmerge, 1, MPI_INTEGER, MPI_SUM, comm_world, mpi_err)

            if (nmerge == 0) then
                call stop_timer(merge_nearest_timer)
                return
            endif

            call collect_from_neighbours

            ! allocate arrays
            allocate(isma(0:nmerge))
            allocate(iclo(nmerge))

            isma = 0
            iclo = 0

            ! Find arrays kc1(lijk) & kc2(lijk) which indicate the parcels in grid cell lijk
            ! through n = node(k), for k = kc1(lijk), kc2(lijk):
            kc1(1) = 1
            do lijk = 1, box%ncell-1
                kc1(lijk+1) = kc1(lijk) + nppc(lijk)
            enddo

            kc2 = kc1 - 1
            j = 0
            do n = 1, n_parcels
                lijk = lloca(n)
                k = kc2(lijk) + 1
                node(k) = n
                kc2(lijk) = k

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

            ! SB: do not use temporary (j) index here, so we will be able to parallelise.
            ! Rather, stop if no nearest parcel found  in surrounding grid boxes
            do m = 1, nmerge
                is = isma(m)
                x_small = parcels%position(1, is)
                y_small = parcels%position(2, is)
                z_small = parcels%position(3, is)
                ! Parcel "is" is small and should be merged; find closest other:
                ix0 = nint(dxi(1) * (x_small - lower(1))) ! ranges from 0 to nx
                iy0 = nint(dxi(2) * (y_small - lower(2))) ! ranges from 0 to ny
                iz0 = nint(dxi(3) * (z_small - lower(3))) ! ranges from 0 to nz

                ! Grid point (ix0, iy0, iz0) is closest to parcel "is"

                dsqmin = product(extent)
                ic = 0

                ! Loop over 8 cells surrounding (ix0, iy0, iz0):
                do iz = max(0, iz0-1), min(box%nz-1, iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do iy = iy0-1, iy0
                        do ix = ix0-1, ix0
                            ! Cell index:
                            gijk = 1 + mod(nx + ix, nx) + nx * mod(ny + iy, ny) + nx * ny * iz
                            ! Search small parcels for closest other:
                            do k = kc1(gijk), kc2(gijk)
                                n = node(k)
                                if (n .ne. is) then
                                    delz = parcels%position(3, n) - z_small
                                    if (delz*delz < dsqmin) then
                                        ! works across periodic edge
                                        delx = get_delx(parcels%position(1, n), x_small)
                                        dely = get_dely(parcels%position(2, n), y_small)
                                        ! Minimise dsqmin
                                        dsq = delx ** 2 + dely ** 2 + delz ** 2
                                        if (dsq < dsqmin) then
                                            dsqmin = dsq
                                            ic = n
                                        endif
                                    endif
                                endif
                            enddo
                        enddo
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

            call stop_timer(merge_nearest_timer)
            call start_timer(merge_tree_resolve_timer)

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

            call stop_timer(merge_tree_resolve_timer)

        end subroutine find_nearest

        subroutine collect_from_neighbours

            ! Identify if parcel in boundary region (which is the halo region of neighbours)
            do iz = box%lo(3), box%hi(3)
                do ix = box%lo(1), box%hi(1)
                    ! north
                    lijk = 1 + ix + box%nx * box%hi(2) + box%nx * box%ny * iz
                    n_halo_small(1) = n_halo_small(1) + nppc(lijk)

                    ! south
                    lijk = 1 + ix + box%nx * box%lo(2) + box%nx * box%ny * iz
                    n_halo_small(2) = n_halo_small(2) + nppc(lijk)
                enddo

                do iy = box%lo(2), box%hi(2)
                    ! west
                    lijk = 1 + box%lo(1) + box%nx * iy + box%nx * box%ny * iz
                    n_halo_small(3) = n_halo_small(3) + nppc(lijk)

                    ! east
                    lijk = 1 + box%hi(1) + box%nx * iy + box%nx * box%ny * iz
                    n_halo_small(4) = n_halo_small(4) + nppc(lijk)
                enddo

                ! north-west
                lijk = 1 + box%lo(1) + box%nx * box%hi(2) + box%nx * box%ny * iz
                n_halo_small(5) = n_halo_small(5) + nppc(lijk)

                ! north-east
                lijk = 1 + box%hi(1) + box%nx * box%hi(2) + box%nx * box%ny * iz
                n_halo_small(6) = n_halo_small(6) + nppc(lijk)

                ! south-west
                lijk = 1 + box%lo(1) + box%nx * box%lo(2) + box%nx * box%ny * iz
                n_halo_small(7) = n_halo_small(7) + nppc(lijk)

                ! south-east
                lijk = 1 + box%hi(1) + box%nx * box%lo(2) + box%nx * box%ny * iz
                n_halo_small(8) = n_halo_small(8) + nppc(lijk)
            enddo

        end subroutine collect_from_neighbours

end module parcel_nearest
