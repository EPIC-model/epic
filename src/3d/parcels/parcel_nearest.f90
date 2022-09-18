!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : pi, f12
    use parcel_container, only : parcels, n_parcels, get_delx, get_dely
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use mpi_communicator
    use mpi_layout
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
    integer :: ic, is, rc, k, m, j, n
    integer :: lijk ! local ijk index
    integer :: gijk ! global ijk index
    integer :: ix, iy, iz, ix0, iy0, iz0
!     integer :: n_halo_small(8)  ! number of small parcels in halo regions

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
        ! @param[out] nmerge the array size of isma and iclo
        ! @post
        !   - isma must be sorted in ascending order
        !   - isma and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(isma, iclo, nmerge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer,              intent(out) :: nmerge
            integer, allocatable              :: rclo(:)

            call start_timer(merge_nearest_timer)

            call nearest_allocate

            nmerge = 0

            !---------------------------------------------------------------------
            ! Initialise search:
            nppc = 0 !nppc(lijk) will contain the number of parcels in grid cell gijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels
                ix = mod(int(dxi(1) * (parcels%position(1, n) - lower(1))), nx)
                iy = mod(int(dxi(2) * (parcels%position(2, n) - lower(2))), ny)
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
            allocate(rclo(nmerge))

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
                l_merged(is)=.false.
                l_merged(ic)=.false.
            enddo

#ifndef NDEBUG
            write(*,*) 'start merging, nmerge='
            write(*,*) nmerge
#endif

            call stop_timer(merge_nearest_timer)

            call resolve_tree(isma, iclo, rclo, nmerge)

        end subroutine find_nearest


        subroutine resolve_tree(isma, iclo, rclo, nmerge)
            integer, intent(inout)         :: isma(:)
            integer, intent(inout)         :: iclo(:)
            integer, intent(inout)         :: rclo(:)
            integer, intent(inout)         :: nmerge
            logical                        :: l_helper
            integer(KIND=MPI_ADDRESS_KIND) :: icm

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

                do m = 1, nmerge
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
                do m = 1, nmerge
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
                do m = 1, nmerge
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
                do m = 1, nmerge
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
            do m = 1, nmerge
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
            do m = 1, nmerge
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


        subroutine collect_from_neighbours
            integer :: n_sends(8)
            integer :: n_recvs(8)

            ! Identify if parcel in boundary region (which is the halo region of neighbours)
            do iz = box%lo(3), box%hi(3)
                do ix = box%lo(1), box%hi(1)
                    ! north
                    lijk = 1 + ix + box%nx * box%hi(2) + box%nx * box%ny * iz
                    n_sends(1) = n_sends(1) + nppc(lijk)

                    ! south
                    lijk = 1 + ix + box%nx * box%lo(2) + box%nx * box%ny * iz
                    n_sends(2) = n_sends(2) + nppc(lijk)
                enddo

                do iy = box%lo(2), box%hi(2)
                    ! west
                    lijk = 1 + box%lo(1) + box%nx * iy + box%nx * box%ny * iz
                    n_sends(3) = n_sends(3) + nppc(lijk)

                    ! east
                    lijk = 1 + box%hi(1) + box%nx * iy + box%nx * box%ny * iz
                    n_sends(4) = n_sends(4) + nppc(lijk)
                enddo

                ! north-west
                lijk = 1 + box%lo(1) + box%nx * box%hi(2) + box%nx * box%ny * iz
                n_sends(5) = n_sends(5) + nppc(lijk)

                ! north-east
                lijk = 1 + box%hi(1) + box%nx * box%hi(2) + box%nx * box%ny * iz
                n_sends(6) = n_sends(6) + nppc(lijk)

                ! south-west
                lijk = 1 + box%lo(1) + box%nx * box%lo(2) + box%nx * box%ny * iz
                n_sends(7) = n_sends(7) + nppc(lijk)

                ! south-east
                lijk = 1 + box%hi(1) + box%nx * box%lo(2) + box%nx * box%ny * iz
                n_sends(8) = n_sends(8) + nppc(lijk)
            enddo

            ! tell your neighbours the number of small parcels in the halo
            call send_to_neighbour(n_sends(NB_NORTH), n_recvs(NB_SOUTH), &
                                   neighbour%north, neighbour%south, NORTH_TAG)

            call send_to_neighbour(n_sends(NB_SOUTH), n_recvs(NB_NORTH), &
                                  neighbour%south, neighbour%north, SOUTH_TAG)

            call send_to_neighbour(n_sends(NB_WEST), n_recvs(NB_EAST), &
                                   neighbour%west, neighbour%east, WEST_TAG)

            call send_to_neighbour(n_sends(NB_EAST), n_recvs(NB_WEST), &
                                  neighbour%east, neighbour%west, EAST_TAG)

            call send_to_neighbour(n_sends(NB_NORTHWEST), n_recvs(NB_SOUTHEAST), &
                                   neighbour%northwest, neighbour%southeast, NORTHWEST_TAG)

            call send_to_neighbour(n_sends(NB_SOUTHEAST), n_recvs(NB_NORTHWEST), &
                                   neighbour%southeast, neighbour%northwest, SOUTHEAST_TAG)

            call send_to_neighbour(n_sends(NB_NORTHEAST), n_recvs(NB_SOUTHWEST), &
                                   neighbour%northeast, neighbour%southwest, NORTHEAST_TAG)

            call send_to_neighbour(n_sends(NB_SOUTHWEST), n_recvs(NB_NORTHEAST), &
                                   neighbour%southwest, neighbour%northeast, SOUTHWEST_TAG)

        end subroutine collect_from_neighbours

        subroutine send_to_neighbour(sendbuf, recvbuf, dest, source, tag)
            integer,          intent(in)  :: sendbuf
            integer,          intent(out) :: recvbuf
            integer,          intent(in)  :: dest, source, tag
            type(MPI_Request)             :: request
            call MPI_Isend(sendbuf, 1, MPI_INT, dest, tag, comm_cart, request, mpi_err)
            call MPI_Request_free(request)
            call MPI_Recv(recvbuf, 1, MPI_INT, source, tag, comm_cart, MPI_STATUS_IGNORE, mpi_err)
        end subroutine send_to_neighbour

end module parcel_nearest
