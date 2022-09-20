!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use dimensions, only : n_dim
    use constants, only : pi, f12
    use parcel_container, only : parcels, n_parcels, get_delx, get_dely
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use mpi_communicator
    use mpi_layout
    use parcel_mpi
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
    logical, allocatable :: l_merged(:) ! indicates parcels merged in first stage

#ifndef NDEBUG
    ! Logicals that are only needed for sanity checks
    logical, allocatable :: l_is_merged(:) ! SANITY CHECK ONLY
    logical, allocatable :: l_small(:)     ! SANITY CHECK ONLY
    logical, allocatable :: l_close(:)     ! SANITY CHECK ONLY
#endif

    logical :: l_continue_iteration, l_do_merge

    !Other variables:
    double precision:: delx, dely, delz, dsq, dsqmin, x_small, y_small, z_small
    integer :: ic, is, rc, k, m, j
    integer :: lijk ! local ijk index
    integer :: gijk ! global ijk index
!     integer :: n_halo_small(8)  ! number of small parcels in halo regions

    integer :: n_small_to_neighbours(8)
    integer :: n_small_neighbours

    type(MPI_Win) :: win_merged, win_avail, win_leaf

#ifndef NDEBUG
    type(MPI_Win) :: win_is_merged, win_small, win_close
#endif

    public :: find_nearest, merge_nearest_timer, merge_tree_resolve_timer

    contains

        subroutine nearest_allocate
            integer, parameter :: disp = 1
            integer (KIND=MPI_ADDRESS_KIND) :: length

            if (.not. allocated(nppc)) then
                allocate(nppc(box%ncell))
                allocate(kc1(box%ncell))
                allocate(kc2(box%ncell))
                allocate(lloca(max_num_parcels))
                allocate(gloca(max_num_parcels))
                allocate(node(max_num_parcels))
                allocate(l_leaf(max_num_parcels))
                allocate(l_available(max_num_parcels))
                allocate(l_merged(max_num_parcels))

                length = sizeof(l_do_merge) * max_num_parcels
                call MPI_Win_create(l_leaf, length, 1, MPI_INFO_NULL, comm_world, win_leaf, mpi_err)
                call MPI_Win_create(l_available, length, 1, MPI_INFO_NULL, comm_world, win_avail, mpi_err)
                call MPI_Win_create(l_merged, length, 1, MPI_INFO_NULL, comm_world, win_merged, mpi_err)

#ifndef NDEBUG
                allocate(l_is_merged(max_num_parcels))
                allocate(l_small(max_num_parcels))
                allocate(l_close(max_num_parcels))

                !     MPI_Win_create(base, size, disp_unit, info, comm, win, ierror)
                !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: base
                !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
                !     INTEGER, INTENT(IN) :: disp_unit
                !     TYPE(MPI_Info), INTENT(IN) :: info
                !     TYPE(MPI_Comm), INTENT(IN) :: comm
                !     TYPE(MPI_Win), INTENT(OUT) :: win
                !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                call MPI_Win_create(l_is_merged, length, 1, MPI_INFO_NULL, comm_world, win_is_merged, mpi_err)
                call MPI_Win_create(l_small, length, 1, MPI_INFO_NULL, comm_world, win_small, mpi_err)
                call MPI_Win_create(l_close, length, 1, MPI_INFO_NULL, comm_world, win_close, mpi_err)
#endif
            endif
        end subroutine nearest_allocate

        ! @param[out] isma indices of small parcels
        ! @param[out] iclo indices of close parcels
        ! @param[out] n_local_merge number of merges
        ! @post
        !   - isma must be sorted in ascending order
        !   - isma and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(isma, iclo, n_local_merge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer,              intent(out) :: n_local_merge
            integer, allocatable              :: rclo(:)
            integer                           :: n_global_merge
            integer                           :: ix, iy, iz, ix0, iy0, iz0
            integer                           :: n

            call start_timer(merge_nearest_timer)

            call nearest_allocate

            n_local_merge = 0
            n_global_merge = 0
            n_small_neighbours = 0 ! number of received small parcels from neighbours
            n_small_to_neighbours(:) = 0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(lijk) will contain the number of parcels in grid cell gijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels
                ix = int(dxi(1) * (parcels%position(1, n) - lower(1)))
                iy = int(dxi(2) * (parcels%position(2, n) - lower(2)))
                iz = min(int(dxi(3) * (parcels%position(3, n) - lower(3))), nz-1)

                ! Cell index of parcel:
                gijk = 1 + ix +     nx * iy +     nx *     ny * iz
                lijk = 1 + ix + box%nx * iy + box%nx * box%ny * iz !This runs from 1 to box%ncell

                ! Accumulate number of parcels in this grid cell:
                nppc(lijk) = nppc(lijk) + 1

                ! Store grid cell that this parcel is in:
                lloca(n) = lijk
                gloca(n) = gijk

                if (parcels%volume(n) < vmin) then
                    n_local_merge = n_local_merge + 1

                    ! If in halo region add to proper location in
                    ! n_small_to_neighbours array. Also, fill send
                    ! buffers with parcel index.
                    call locate_small_parcel(n, ix, iy)
                endif
            enddo

            call MPI_Allreduce(n_local_merge, n_global_merge, 1, MPI_INTEGER, MPI_SUM, comm_world, mpi_err)

            if (n_global_merge == 0) then
                call stop_timer(merge_nearest_timer)
                return
            endif

            ! send small parcels to neighbours
            call send_small_to_neighbours(north_buf, neighbour%north, neighbour%south, &
                                          NORTH_TAG, NB_NORTH, NB_SOUTH)

            call send_small_to_neighbours(east_buf, neighbour%east, neighbour%west, &
                                          EAST_TAG, NB_EAST, NB_WEST)

            call send_small_to_neighbours(south_buf, neighbour%south, neighbour%north, &
                                          SOUTH_TAG, NB_SOUTH, NB_NORTH)

            call send_small_to_neighbours(west_buf, neighbour%west, neighbour%east, &
                                          WEST_TAG, NB_WEST, NB_EAST)

            call send_small_to_neighbours(northwest_buf, neighbour%northwest, neighbour%northeast, &
                                          NORTHWEST_TAG, NB_NORTHWEST, NB_NORTHEAST)

            call send_small_to_neighbours(northeast_buf, neighbour%northeast, neighbour%northwest, &
                                          NORTHEAST_TAG, NB_NORTHEAST, NB_NORTHWEST)

            call send_small_to_neighbours(southwest_buf, neighbour%southwest, neighbour%southeast, &
                                          SOUTHWEST_TAG, NB_SOUTHWEST, NB_SOUTHEAST)

            call send_small_to_neighbours(southeast_buf, neighbour%southeast, neighbour%southwest, &
                                          SOUTHEAST_TAG, NB_SOUTHEAST, NB_SOUTHWEST)


            if ((n_local_merge == 0) .and. (n_small_neighbours > 0)) then
                ! This rank has no local small parcels. To be not involved in any merging,
                ! it must send all its close parcels to its neighbours owning the small
                ! parcel. Although this rank is not involved in any merging, it must nevertheless
                ! call the tree resolving routine as there is global synchronisation taking place.

!                 call send_close_to_neighbours

                call resolve_tree(isma, iclo, rclo, n_local_merge)

                return

            else if (n_local_merge + n_small_neighbours == 0) then
                ! This rank is involved in no merging process. Nevertheless, it must
                ! call the tree resolving routine as there is global synchronisation taking place
                call resolve_tree(isma, iclo, rclo, n_local_merge)

                return
            endif

            ! allocate arrays
            allocate(isma(0:n_local_merge + n_small_neighbours))
            allocate(iclo(n_local_merge + n_small_neighbours))
            allocate(rclo(n_local_merge + n_small_neighbours))


            isma = 0
            iclo = 0
            rclo = -1

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
                l_is_merged(n)=.false.! SANITY CHECK ONLY
                l_small(n)=.false. ! SANITY CHECK ONLY
                l_close(n)=.false. ! SANITY CHECK ONLY
#endif
            enddo

            !---------------------------------------------------------------------
            ! Now find the nearest grid point to each small parcel (to be merged)
            ! and search over the surrounding 8 grid cells for the closest parcel:

            ! SB: do not use temporary (j) index here, so we will be able to parallelise.
            ! Rather, stop if no nearest parcel found  in surrounding grid boxes
            do m = 1, n_local_merge
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
                l_merged(is)=.false.
                l_merged(ic)=.false.
            enddo

#ifndef NDEBUG
            write(*,*) 'start merging, nmerge='
            write(*,*) n_local_merge
#endif

            call stop_timer(merge_nearest_timer)

            call resolve_tree(isma, iclo, rclo, n_local_merge)

        end subroutine find_nearest


        subroutine resolve_tree(isma, iclo, rclo, n_local_merge)
            integer, intent(inout)         :: isma(:)
            integer, intent(inout)         :: iclo(:)
            integer, intent(inout)         :: rclo(:)
            integer, intent(inout)         :: n_local_merge
            logical                        :: l_helper
            integer(KIND=MPI_ADDRESS_KIND) :: icm
#ifndef NDEBUG
            integer                        :: n
#endif

            call start_timer(merge_tree_resolve_timer)

            ! First implementation of iterative algorithm

            ! Could be improved by keeping track of "finalised mergers"
            ! But going for more readable implementation here first.

            ! First, iterative, stage
            l_continue_iteration = .true.

            do while(l_continue_iteration)
                l_continue_iteration = .false.
                ! reset relevant properties for candidate mergers

                call MPI_Win_fence(0, win_avail, mpi_err)

                do m = 1, n_local_merge
                    is = isma(m)
                    ! only consider links that still may be merging
                    ! reset relevant properties
                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)
                        l_leaf(is) = .true.

                        if (rc == mpi_rank) then
                            l_available(ic) = .true.
                        else
                            !     MPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                            !         target_disp, target_count, target_datatype, win, ierror)
                            !     TYPE(*), DIMENSION(..), INTENT(IN), ASYNCHRONOUS :: origin_addr
                            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
                            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
                            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
                            !     TYPE(MPI_Win), INTENT(IN) :: win
                            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                            l_helper = .true.
                            icm = ic
                            call MPI_Put(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                         MPI_LOGICAL, win_avail, mpi_err)
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_avail, mpi_err)

                call MPI_Win_fence(0, win_leaf, mpi_err)

                ! determine leaf parcels
                do m = 1, n_local_merge
                    is = isma(m)
                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        if (rc == mpi_rank) then
                            l_leaf(ic) = .false.
                        else
                            l_helper = .false.
                            icm = ic
                            call MPI_Put(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                         MPI_LOGICAL, win_leaf, mpi_err)
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_leaf, mpi_err)

                call MPI_Win_fence(0, win_avail, mpi_err)

                ! filter out parcels that are "unavailable" for merging
                do m = 1, n_local_merge
                    is = isma(m)
                    if (.not. l_merged(is)) then
                        if (.not. l_leaf(is)) then
                            ic = iclo(m)
                            rc = rclo(m)

                            if (rc == mpi_rank) then
                                l_available(ic) = .false.
                            else
                                l_helper = .false.
                                icm = ic
                                call MPI_Put(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                             MPI_LOGICAL, win_avail, mpi_err)
                            endif
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_avail, mpi_err)

                call MPI_Win_fence(0, win_merged, mpi_err)

                ! identify mergers in this iteration
                do m = 1, n_local_merge
                    is = isma(m)
                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        l_helper = .false.
                        if (rc == mpi_rank) then
                            l_helper = l_available(ic)
                        else
                            !     MPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                            !         target_disp, target_count, target_datatype, win, ierror)
                            !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: origin_addr
                            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
                            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
                            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
                            !     TYPE(MPI_Win), INTENT(IN) :: win
                            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                            icm = ic
                            call MPI_Get(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                         MPI_LOGICAL, win_avail, mpi_err)
                        endif

                        if (l_leaf(is) .and. l_helper) then
                            l_continue_iteration = .true. ! merger means continue iteration
                            l_merged(is) = .true.

                            if (rc == mpi_rank) then
                                l_merged(ic) = .true.
                            else
                                l_helper = .true.
                                icm = ic
                                call MPI_Put(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                             MPI_LOGICAL, win_merged, mpi_err)
                            endif
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_merged, mpi_err)

                ! Performance improvement: We actually only need to synchronize with neighbours
                call MPI_Allreduce(MPI_IN_PLACE, l_continue_iteration, 1, MPI_LOGICAL, &
                                   MPI_LOR, comm_world, mpi_err)
            enddo


            call MPI_Win_fence(0, win_avail, mpi_err)

            ! Second stage, related to dual links
            do m = 1, n_local_merge
                is = isma(m)
                if (.not. l_merged(is)) then
                    if (l_leaf(is)) then ! set in last iteration of stage 1
                        ic = iclo(m)
                        rc = rclo(m)

                        if (rc == mpi_rank) then
                            l_available(ic) = .true.
                        else
                            l_helper = .true.
                            icm = ic
                            call MPI_Put(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                         MPI_LOGICAL, win_avail, mpi_err)
                        endif
                    endif
                endif
            enddo

            call MPI_Win_fence(0, win_avail, mpi_err)

            call MPI_Win_fence(0, win_avail, mpi_err)

            ! Second stage (hard to parallelise with openmp?)
            j = 0
            do m = 1, n_local_merge
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                l_do_merge = .false.
                if (l_merged(is) .and. l_leaf(is)) then
                    ! previously identified mergers: keep
                    l_do_merge = .true.
#ifndef NDEBUG
                    ! sanity check on first stage mergers
                    ! parcel cannot be both initiator and receiver in stage 1
                    if (l_leaf(ic)) then
                        write(*,*) 'first stage error'
                    endif
#endif
                elseif (.not. l_merged(is)) then
                    if (l_leaf(is)) then
                        ! links from leafs
                        l_do_merge = .true.
                    elseif (.not. l_available(is)) then
                        ! Above means parcels that have been made 'available' do not keep outgoing links

                        if (rc == mpi_rank) then
                            l_helper = l_available(ic)
                        else
                            call MPI_Get(l_helper, 1, MPI_LOGICAL, rc, icm, max_num_parcels, &
                                         MPI_LOGICAL, win_avail, mpi_err)
                        endif

                        if (l_helper) then
                            ! merge this parcel into ic along with the leaf parcels
                            l_do_merge = .true.
                        else
                            ! isolated dual link
                            ! Don't keep current link
                            ! But make small parcel available so other parcel can merge with it
                            ! THIS NEEDS THINKING ABOUT A PARALLEL IMPLEMENTATION
                            ! This could be based on the other parcel being outside the domain
                            ! And a "processor order"
                            l_available(is) = .true.
                        endif
                    endif
              endif

              call MPI_Win_fence(0, win_avail, mpi_err)

              if (l_do_merge) then
                   j = j + 1
                   isma(j) = is
                   iclo(j) = ic
                   rclo(j) = rc
#ifndef NDEBUG
                   l_is_merged(is)=.true.
                   l_is_merged(ic)=.true.
                   l_small(is)=.true.
                   l_close(ic)=.true.
#endif
              end if
            end do
            n_local_merge = j

#ifndef NDEBUG
            write(*,*) 'after second stage, nmerge='
            write(*,*) n_local_merge
            write(*,*) 'finished'

            ! MORE SANITY CHECKS
            ! CHECK ISMA ORDER
            do m = 1, n_local_merge
              if(.not. (isma(m)>isma(m-1))) then
                write(*,*) 'isma order broken'
              end if
            end do

            ! 1. CHECK RESULTING MERGERS
            do m = 1, n_local_merge
              if(.not. l_is_merged(isma(m))) write(*,*) 'merge_error: isma(m) not merged, m=', m
              if(.not. l_is_merged(iclo(m))) write(*,*) 'merge_error: iclo(m) not merged, m=', m
              if(.not. l_small(isma(m))) write(*,*) 'merge_error: isma(m) not marked as small, m=', m
              if(.not. l_close(iclo(m))) write(*,*) 'merge_error: iclo(m) not marked as close, m=', m
              if(l_close(isma(m))) write(*,*) 'merge_error: isma(m) both small and close, m=', m
              if(l_small(iclo(m))) write(*,*) 'merge_error: iclo(m) both small and close, m=', m
            end do

            ! 2. CHECK MERGING PARCELS
            do n = 1, n_parcels
              if (parcels%volume(n) < vmin) then
                if(.not. l_is_merged(n)) write(*,*) 'merge_error: parcel n not merged (should be), n=', n
                if(.not. (l_small(n) .or. l_close(n))) write(*,*) 'merge_error: parcel n not small or close (should be), n=', n
                if(l_small(n) .and. l_close(n)) write(*,*) 'merge_error: parcel n both small and close, n=', n
              else
                if(l_small(n)) write(*,*) 'merge_error: parcel n small (should not be), n=', n
                if(l_is_merged(n) .and. (.not. l_close(n))) write(*,*) 'merge_error: parcel n merged (should not be), n=', n
              end if
            enddo
#endif

            call stop_timer(merge_tree_resolve_timer)

        end subroutine resolve_tree

        subroutine locate_small_parcel(n, ix, iy)
            integer, intent(in) :: n, ix, iy
            integer             :: m

            ! check lower x-direction
            if (ix == box%lo(1)) then

                if (iy == box%lo(2)) then     ! check southwest corner
                    ! neighbours: west, south and southwest
                    m = n_small_to_neighbours(NB_SOUTH) + 1
                    south_buf(m) = n
                    n_small_to_neighbours(NB_SOUTH) = m

                    m = n_small_to_neighbours(NB_SOUTHWEST) + 1
                    southwest_buf(m) = n
                    n_small_to_neighbours(NB_SOUTHWEST) = m

                else if (iy == box%hi(2)) then
                    ! neighbours: west, north and northwest
                    m = n_small_to_neighbours(NB_NORTH) + 1
                    north_buf(m) = n
                    n_small_to_neighbours(NB_NORTH) = m

                    m = n_small_to_neighbours(NB_NORTHWEST) + 1
                    northwest_buf(m) = n
                    n_small_to_neighbours(NB_NORTHWEST) = m
                endif

                ! neighbour: west
                m = n_small_to_neighbours(NB_WEST) + 1
                west_buf(m) = n
                n_small_to_neighbours(NB_WEST) = m
            endif

            ! check upper x-direction (use >= as a parcel might be directly on the cell edge)
            if (ix >= box%hi(1)) then

                if (iy == box%lo(2)) then     ! check southeast corner
                    ! neighbours: east, south and southeast
                    m = n_small_to_neighbours(NB_SOUTH) + 1
                    south_buf(m) = n
                    n_small_to_neighbours(NB_SOUTH) = m

                    m = n_small_to_neighbours(NB_SOUTHEAST) + 1
                    southeast_buf(m) = n
                    n_small_to_neighbours(NB_SOUTHEAST) = m

                else if (iy == box%hi(2)) then
                    ! neighbours: east, north and northeast
                    m = n_small_to_neighbours(NB_NORTH) + 1
                    north_buf(m) = n
                    n_small_to_neighbours(NB_NORTH) = m

                    m = n_small_to_neighbours(NB_NORTHEAST) + 1
                    northeast_buf(m) = n
                    n_small_to_neighbours(NB_NORTHEAST) = m
                endif

                ! neighbours: east
                m = n_small_to_neighbours(NB_EAST) + 1
                east_buf(m) = n
                n_small_to_neighbours(NB_EAST) = m
            endif

        end subroutine locate_small_parcel

        ! Send position attribute of small parcels in boundary region
        ! to neighbours. We only need to send the position as this is the
        ! only attribute needed to figure out with whom a parcel might merge.
        subroutine send_small_to_neighbours(pid, dest, source, tag, sendloc, recvloc)
            integer,                      intent(in) :: pid(:)
            integer, intent(in)                      :: dest, source, tag
            integer, intent(in)                      :: sendloc, recvloc
            double precision, allocatable            :: sendbuf(:), recvbuf(:)
            integer                                  :: send_size, recv_size
            integer                                  :: n, i, j, sendcount, recvcount, k
            type(MPI_Request)                        :: request

            sendcount = n_small_to_neighbours(sendloc)
            recvcount = n_small_to_neighbours(recvloc)

            send_size = n_dim * n_small_to_neighbours(sendloc)
            recv_size = n_dim * n_small_to_neighbours(recvloc)

            allocate(sendbuf(send_size))
            allocate(recvbuf(recv_size))

            ! pack parcel positions to send buffer
            do n = 1, sendcount
                i = 1 + (n-1) * n_dim
                j = n * n_dim
                sendbuf(i:j) = parcels%position(:, pid(n))
            enddo

            call MPI_Isend(sendbuf, send_size, MPI_DOUBLE, dest, tag, comm_cart, request, mpi_err)
            call MPI_Request_free(request)

            call MPI_Recv(recvbuf, recv_size, MPI_DOUBLE, source, tag, &
                          comm_cart, MPI_STATUS_IGNORE, mpi_err)


            ! unpack parcel positions to recv buffer
            do n = 1, recvcount
                i = 1 + (n-1) * n_dim
                j = n * n_dim
                k = n_parcels + n
                parcels%position(:, k) = recvbuf(i:j)
            enddo

            n_small_neighbours = n_small_neighbours + recvcount

            deallocate(sendbuf)
            deallocate(recvbuf)

        end subroutine send_small_to_neighbours

end module parcel_nearest
