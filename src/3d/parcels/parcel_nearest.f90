!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : zero !pi, f12
    use parcel_container, only : parcels, n_parcels, get_delx, get_dely
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use timer, only : start_timer, stop_timer
    use mpi_communicator
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    integer:: merge_nearest_timer, merge_tree_resolve_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:), kc2(:)
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

    !FIXME
    integer :: n_parcel_sends(8)
    integer :: southwest_pid(1), southeast_pid(1), northwest_pid(1), northeast_pid(1)
    integer :: east_pid(1), west_pid(1), south_pid(1), north_pid(1)

    public :: find_nearest, merge_nearest_timer, merge_tree_resolve_timer

    type :: t_pair
        integer :: rank     ! MPI rank
        integer :: pid      ! parcel id
    end type t_pair

    type(t_pair), allocatable :: nsma(:)

    contains

        subroutine nearest_allocate
            if (.not. allocated(nppc)) then
                allocate(nppc(box%halo_ncell))
                allocate(kc1(box%halo_ncell))
                allocate(kc2(box%halo_ncell))
                allocate(loca(max_num_parcels))
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
        end subroutine nearest_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! @param[out] isma indices of small parcels
        ! @param[out] iclo indices of close parcels
        ! @param[out] n_local_merge the array size of isma and iclo
        ! @post
        !   - isma must be sorted in ascending order
        !   - isma and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(isma, iclo, n_local_merge)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer,              intent(out) :: n_local_merge
            integer                           :: n_global_merge
            integer                           :: n_received_small
            integer                           :: ijk, n, ix, iy, k, j
            integer, allocatable              :: rclo(:)

            call start_timer(merge_nearest_timer)

            call nearest_allocate

            !---------------------------------------------------------------------
            ! Initialise search:
            n_local_merge = 0
            n_global_merge = 0
            n_received_small = 0
            nppc = 0 !nppc(ijk) will contain the number of parcels in grid cell ijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels

                call parcel_to_local_cell_index(n, ix, iy)

                if (parcels%volume(n) < vmin) then
                    n_local_merge = n_local_merge + 1

                    ! If a small parcel is in a boundary cell, a duplicate must
                    ! be sent to the neighbour rank. This call checks if the parcel
                    ! must be sent and fills the send buffers.
                    call locate_parcel_in_boundary_cell(n, ix, iy)
                endif
            enddo

            call MPI_Allreduce(n_local_merge, n_global_merge, 1, MPI_INTEGER, MPI_SUM, comm%world, comm%err)

            if (n_global_merge == 0) then
                call stop_timer(merge_nearest_timer)
                return
            endif

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute of small parcels in boundary region
            ! to neighbours. We only need to send the position as this is the
            ! only attribute needed to figure out with whom a parcel might merge.
!             ! Note: n_parcels gets increased with this operation, we must
!             !       revert n_parcels to the previous value as we keep track of
!             !       remote parcels differently.
!             call parcel_send(n_received_small, nsma) !FIXME
!             n_parcels = n_parcels - n_received_small

            ! We must also assign incoming small parcels to cells
            ! Note: We must apply a shift for parcels communicated
            !       across a periodic boundary.
            do n = n_parcels + 1, n_parcels + n_received_small
                ! global cell index in x and y
                ix = int(dxi(1) * (parcels%position(1, n) - lower(1)))
                iy = int(dxi(2) * (parcels%position(2, n) - lower(2)))

                ! Check if parcel is outside *this* subdomain (including halo region) and
                ! shift if necessary. Across a periodic boundary we must shift a parcel
                ! into the halo region of *this* subdomain.

                ! If (/ix, iy/) < box%hlo(1:2) is true, then add global domain
                ! extent to parcel position.
                parcels%position(1:2, n) = parcels%position(1:2, n) &
                                         + merge(extent(1:2), (/zero, zero/), (/ix, iy/) < box%hlo(1:2))

                ! If (/ix, iy/) > box%hhi(1:2) is true, then subtract global domain
                ! extent from parcel position.
                parcels%position(1:2, n) = parcels%position(1:2, n) &
                                         - merge(extent(1:2), (/zero, zero/), (/ix, iy/) > box%hhi(1:2))

                ! Now we can safely assign the local cell index:
                call parcel_to_local_cell_index(n, ix, iy)
            enddo

            ! There are 4 cases:
            !   - n_local_merge = 0 and n_received_small = 0 --> This rank has no small parcels.
            !   - n_local_merge > 0 and n_received_small = 0 --> This rank has only local small parcels.
            !   - n_local_merge = 0 and n_received_small > 0 --> This rank has only remote small parcels.
            !   - n_local_merge > 0 and n_received_small > 0 --> This rank has local and remote small parcels.
            if (n_local_merge + n_received_small == 0) then
                ! Although *this* rank has no small parcels and is therefore not involved
                ! in any merging, it must nonetheless call the tree resolving routine as
                ! there is global synchronisation taking place.
                call resolve_tree(isma, iclo, rclo, n_local_merge)

                return
            endif

            ! allocate arrays
            allocate(isma(0:n_local_merge + n_received_small))
            allocate(iclo(n_local_merge + n_received_small))
            allocate(rclo(n_local_merge + n_received_small))

            isma = 0
            iclo = 0
            rclo = -1

            ! Find arrays kc1(ijk) & kc2(ijk) which indicate the parcels in grid cell ijk
            ! through n = node(k), for k = kc1(ijk), kc2(ijk):
            kc1(1) = 1
            do ijk = 1, box%halo_ncell-1
                kc1(ijk+1) = kc1(ijk) + nppc(ijk)
            enddo

            kc2 = kc1 - 1
            j = 0
            do n = 1, n_parcels
                ijk = loca(n)
                k = kc2(ijk) + 1
                node(k) = n
                kc2(ijk) = k

                if (parcels%volume(n) < vmin) then
                    j = j + 1
                    isma(j) = n
                endif
#ifndef NDEBUG
                l_merged(n) = .false.! SANITY CHECK ONLY
                l_small(n) = .false. ! SANITY CHECK ONLY
                l_close(n) = .false. ! SANITY CHECK ONLY
#endif
            enddo

            !---------------------------------------------------------------------
            ! Determine locally closest parcel:

            call find_closest_parcel_locally(isma, iclo, rclo, n_local_merge)

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:

            call find_closest_parcel_globally(isma, iclo, rclo, n_local_merge, n_received_small)

#ifndef NDEBUG
            write(*,*) 'start merging, n_local_merge='
            write(*,*) n_local_merge
#endif

            call stop_timer(merge_nearest_timer)

            call resolve_tree(isma, iclo, rclo, n_local_merge)

        end subroutine find_nearest

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !!@pre Assumes a parcel is in the local domain including halo cells
        !      (in x and y).
        subroutine parcel_to_local_cell_index(n, ix, iy)
            integer, intent(in)  :: n
            integer, intent(out) :: ix, iy
            integer              :: iz, ijk

            ix =     int(dxi(1) * (parcels%position(1, n) - box%halo_lower(1)))
            iy =     int(dxi(2) * (parcels%position(2, n) - box%halo_lower(2)))
            iz = min(int(dxi(3) * (parcels%position(3, n) - box%lower(3))), nz-1)

            ! Cell index of parcel:
            !   This runs from 1 to halo_ncell where
            !   halo_ncell includes halo cells
            ijk = 1 + ix + box%halo_size(1) * iy + box%halo_size(1) * box%halo_size(2) * iz

            ! Accumulate number of parcels in this grid cell:
            nppc(ijk) = nppc(ijk) + 1

            ! Store grid cell that this parcel is in:
            loca(n) = ijk

        end subroutine parcel_to_local_cell_index

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Find the nearest grid point to each small parcel (to be merged)
        ! and search over the surrounding 8 grid cells for the closest parcel.
        ! This operation is performed locally.
        subroutine find_closest_parcel_locally(isma, iclo, rclo, n_local_merge)
            integer, intent(inout) :: isma(0:)
            integer, intent(inout) :: iclo(:)
            integer, intent(inout) :: rclo(:)
            integer, intent(in)    :: n_local_merge
            integer                :: ix0, iy0, iz0, ijk, n, ix, iy, iz, m, k, ic, is
            double precision       :: delx, dely, delz, dsq, dsqmin, x_small, y_small, z_small

            ! SB: do not use temporary (j) index here, so we will be able to parallelise.
            ! Rather, stop if no nearest parcel found  in surrounding grid boxes
            do m = 1, n_local_merge
                is = isma(m)
                x_small = parcels%position(1, is)
                y_small = parcels%position(2, is)
                z_small = parcels%position(3, is)
                ! Parcel "is" is small and should be merged; find closest other:
                ix0 = nint(dxi(1) * (x_small - box%halo_lower(1))) ! ranges from 0 to box%halo_size(1)
                iy0 = nint(dxi(2) * (y_small - box%halo_lower(2))) ! ranges from 0 to box%halo_size(2)
                iz0 = nint(dxi(3) * (z_small - box%lower(3)))      ! ranges from 0 to nz

                ! Grid point (ix0, iy0, iz0) is closest to parcel "is"

                dsqmin = product(extent)
                ic = 0

                ! Loop over 8 cells surrounding (ix0, iy0, iz0):
                do iz = max(0, iz0-1), min(nz-1, iz0) !=> iz=0 for iz0=0 & iz=nz-1 for iz0=nz
                    do iy = iy0-1, iy0
                        do ix = ix0-1, ix0
                            ! Cell index:
                            ijk = 1 + ix                                    &
                                + box%halo_size(1) * iy                     &
                                + box%halo_size(1) * box%halo_size(2) * iz

                            ! Search small parcels for closest other:
                            do k = kc1(ijk), kc2(ijk)
                                n = node(k)
                                if (n .ne. is) then
                                    delz = parcels%position(3, n) - z_small
                                    if (delz * delz < dsqmin) then
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
                    call mpi_exit_on_error('Merge error: no near neighbour found.')
                endif

                ! Store the MPI rank and the index of the parcel to be potentially merged with:
                isma(m) = is
                iclo(m) = ic
                rclo(m) = merge(nsma(ic)%rank, comm%rank, ic > n_parcels)
                l_first_merged(is) = .false.
                l_first_merged(ic) = .false.
            enddo
        end subroutine find_closest_parcel_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! For all small parcels in boundary cells we must now determine the
        ! closest parcel globally (that is we must compare with the neighbour rank).
        ! For this purpose we send the information of the small parcels that we received
        ! back to the original ranks which then figure out the closest parcel.
        subroutine find_closest_parcel_globally(isma, iclo, rclo, n_local_merge, n_received_small)
            integer, intent(inout) :: isma(0:)
            integer, intent(inout) :: iclo(:)
            integer, intent(inout) :: rclo(:)
            integer, intent(in)    :: n_local_merge, n_received_small


            integer,          dimension(:), pointer :: pid
            double precision, dimension(:), pointer :: sendbuf
            double precision, allocatable           :: recvbuf(:)
            type(MPI_Request)                       :: request
            type(MPI_Status)                        :: recv_status, send_status
            integer                                 :: recv_size, send_size
            integer                                 :: tag, source, recvcount, n, k!, k, i ,j, l
            integer                                 :: n_sends(8)

            !------------------------------------------------------------------
            ! Figure out how many small parcels we received:
            n_sends = 0
            do k = 1, n_received_small
                n = get_neighbour_from_rank(nsma(n_parcels + k)%rank)
                n_sends(n) = n_sends(n) + 1
            enddo

            !------------------------------------------------------------------
            ! Communicate with neighbours:
            do n = 1, 8

                ! we send the distance and the parcel index
                send_size = n_sends(n) * 2

                allocate(sendbuf(send_size))

                do k = 1, send_size
                    sendbuf(2*k-1) = nsma(
                    sendbuf(2*k)   =
                enddo

                call MPI_Isend(sendbuf,                 &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               request,                 &
                               comm%err)

                call mpi_check_for_error(&
                    "in MPI_Isend of parcel_nearest::find_closest_parcel_globally.")

                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                allocate(recvbuf(recv_size))

                call MPI_Recv(recvbuf,                  &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                call mpi_check_for_error(&
                    "in MPI_Recv of parcel_nearest::find_closest_parcel_globally.")

                call MPI_Wait(request, send_status, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Wait of parcel_nearest::find_closest_parcel_globally.")
            enddo

        end subroutine find_closest_parcel_globally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine resolve_tree(isma, iclo, rclo, n_local_merge)
            integer, intent(inout)         :: isma(0:)
            integer, intent(inout)         :: iclo(:)
            integer, intent(inout)         :: rclo(:)
            integer, intent(inout)         :: n_local_merge
            integer                        :: ic, is, m, j
#ifndef NDEBUG
            integer                        :: n
#endif

            rclo = 0 ! FIXME

            call start_timer(merge_tree_resolve_timer)

            ! First implementation of iterative algorithm

            ! Could be improved by keeping track of "finalised mergers"
            ! But going for more readable implementation here first.

            ! First, iterative, stage
            l_continue_iteration=.true.
            do while(l_continue_iteration)
              l_continue_iteration=.false.
              ! reset relevant properties for candidate mergers
              do m = 1, n_local_merge
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
              do m = 1, n_local_merge
                is = isma(m)
                if(.not. l_first_merged(is)) then
                  ic = iclo(m)
                  l_leaf(ic)=.false.
                end if
              end do
              ! PARALLEL COMMUNICATION WILL BE NEEDED HERE
              ! filter out parcels that are "unavailable" for merging
              do m = 1, n_local_merge
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
              do m = 1, n_local_merge
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
            do m = 1, n_local_merge
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
            do m = 1, n_local_merge
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
            n_local_merge = j

#ifndef NDEBUG
            write(*,*) 'after second stage, n_local_merge='
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
                    if (.not. l_merged(n)) then
                        write(*,*) 'merge_error: parcel n not merged (should be), n=', n
                    endif
                    if (.not. (l_small(n) .or. l_close(n))) then
                        write(*,*) 'merge_error: parcel n not small or close (should be), n=', n
                    endif
                    if(l_small(n) .and. l_close(n)) then
                        write(*,*) 'merge_error: parcel n both small and close, n=', n
                    endif
                else
                    if (l_small(n)) then
                        write(*,*) 'merge_error: parcel n small (should not be), n=', n
                    endif
                    if (l_merged(n) .and. (.not. l_close(n))) then
                        write(*,*) 'merge_error: parcel n merged (should not be), n=', n
                    endif
                endif
            enddo
#endif
            call stop_timer(merge_tree_resolve_timer)
        end subroutine resolve_tree

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine locate_parcel_in_boundary_cell(n, ix, iy)
            integer, intent(in) :: n, ix, iy
            integer             :: m

            ! check lower x-direction
            if (ix == box%lo(1)) then

                if (iy == box%lo(2)) then
                    ! parcel in southwest corner with
                    ! neighbours: west, south and southwest
                    m = n_parcel_sends(MPI_SOUTHWEST) + 1
                    southwest_pid(m) = n
                    n_parcel_sends(MPI_SOUTHWEST) = m

                else if (iy == box%hi(2)) then
                    ! parcel in northwest corner with
                    ! neighbours: west, north and northwest
                    m = n_parcel_sends(MPI_NORTHWEST) + 1
                    northwest_pid(m) = n
                    n_parcel_sends(MPI_NORTHWEST) = m
                endif

                ! neighbour: west
                m = n_parcel_sends(MPI_WEST) + 1
                west_pid(m) = n
                n_parcel_sends(MPI_WEST) = m
            endif

            ! check upper x-direction (use >= as a parcel might be directly on the cell edge)
            if (ix >= box%hi(1)) then

                if (iy == box%lo(2)) then
                    ! parcel in southeast corner with
                    ! neighbours: east, south and southeast
                    m = n_parcel_sends(MPI_SOUTHEAST) + 1
                    southeast_pid(m) = n
                    n_parcel_sends(MPI_SOUTHEAST) = m

                else if (iy == box%hi(2)) then
                    ! parcel in northeast corner with
                    ! neighbours: east, north and northeast
                    m = n_parcel_sends(MPI_NORTHEAST) + 1
                    northeast_pid(m) = n
                    n_parcel_sends(MPI_NORTHEAST) = m
                endif

                ! neighbours: east
                m = n_parcel_sends(MPI_EAST) + 1
                east_pid(m) = n
                n_parcel_sends(MPI_EAST) = m
            endif

            if (iy == box%lo(2)) then
                m = n_parcel_sends(MPI_SOUTH) + 1
                south_pid(m) = n
                n_parcel_sends(MPI_SOUTH) = m
            endif

            if (iy >= box%hi(2)) then
                m = n_parcel_sends(MPI_NORTH) + 1
                north_pid(m) = n
                n_parcel_sends(MPI_NORTH) = m
            endif
        end subroutine locate_parcel_in_boundary_cell

end module parcel_nearest
