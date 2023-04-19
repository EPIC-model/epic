!==============================================================================
!               Finds the parcels nearest to every "small" parcel
!
!   In a first step each MPI rank sends the positions of all parcels in
!   its domain boundary region to the corresponding neighbours. Each rank
!   can then find the local closest parcel to a small parcel. After that
!   the information about the closest parcel is then sent back to the original
!   MPI rank that then determines the closest parcel globally. Once this is
!   done, the tree resolving algorithm is executed where we use one-sided
!   communication with passive target synchronisation. Finally, the small
!   parcels with a closest parcel on another sub-domain must be sent to the
!   appropriate MPI rank.
!==============================================================================
module parcel_nearest
    use constants, only : zero
    use parcel_container, only : parcels            &
                               , n_parcels          &
                               , get_delx           &
                               , get_dely           &
                               , n_par_attrib       &
                               , n_total_parcels    &
                               , parcel_serialize   &
                               , parcel_deserialize &
                               , parcel_delete
    use parameters, only : dx, dxi, vcell, hli, lower, extent       &
                         , ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use timer, only : start_timer, stop_timer
    use mpi_communicator
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_error, mpi_check_for_message
    use iso_c_binding, only : c_sizeof
    use fields, only : get_index_periodic
    use parcel_mpi, only : n_parcel_sends               &
                         , north_pid                    &
                         , south_pid                    &
                         , west_pid                     &
                         , east_pid                     &
                         , northwest_pid                &
                         , northeast_pid                &
                         , southwest_pid                &
                         , southeast_pid                &
                         , north_buf                    &
                         , south_buf                    &
                         , west_buf                     &
                         , east_buf                     &
                         , northwest_buf                &
                         , northeast_buf                &
                         , southwest_buf                &
                         , southeast_buf                &
                         , invalid                      &
                         , allocate_parcel_buffers      &
                         , deallocate_parcel_buffers    &
                         , allocate_parcel_id_buffers   &
                         , deallocate_parcel_id_buffers &
                         , get_parcel_buffer_ptr
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
    logical, allocatable, asynchronous :: l_leaf(:)
    logical, allocatable, asynchronous :: l_available(:)
    logical, allocatable, asynchronous :: l_merged(:)    ! indicates parcels merged in first stage

    integer              :: n_neighbour_small(8)  ! number of small parcels received
    integer              :: small_recv_order(8)   ! receive order of small parcels
    integer              :: n_remote_small        ! sum(n_neighbour_small)
    integer, allocatable :: rsma(:)               ! rank of small parcel received (accessed with *n*)
    integer, allocatable :: pidsma(:)             ! index of remote small parcel (accessed with *n*)
    integer, allocatable :: midsma(:)             ! *m* in *isma(m)*

    type(MPI_Win) :: win_merged, win_avail, win_leaf

    public :: find_nearest, merge_nearest_timer, merge_tree_resolve_timer

    contains

        subroutine nearest_allocate
            integer (KIND=MPI_ADDRESS_KIND) :: win_size
            logical                         :: dummy
            integer                         :: disp_unit = c_sizeof(dummy) ! size of logical in bytes

            if (.not. allocated(nppc)) then
                allocate(nppc(box%halo_ncell))
                allocate(kc1(box%halo_ncell))
                allocate(kc2(box%halo_ncell))
                allocate(loca(max_num_parcels))
                allocate(node(max_num_parcels))
                allocate(l_leaf(max_num_parcels))
                allocate(l_available(max_num_parcels))
                allocate(l_merged(max_num_parcels))

                call MPI_Sizeof(dummy, disp_unit, comm%err)

                ! size of RMA window in bytes
                win_size = disp_unit * max_num_parcels

                !     MPI_Win_create(base, size, disp_unit, info, comm, win, ierror)
                !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: base
                !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: size
                !     INTEGER, INTENT(IN) :: disp_unit
                !     TYPE(MPI_Info), INTENT(IN) :: info
                !     TYPE(MPI_Comm), INTENT(IN) :: comm
                !     TYPE(MPI_Win), INTENT(OUT) :: win
                !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                call MPI_Win_create(l_leaf(1:max_num_parcels),  &
                                    win_size,                   &
                                    disp_unit,                  &
                                    MPI_INFO_NULL,              &
                                    comm%world,                 &
                                    win_leaf,                   &
                                    comm%err)
                call MPI_Win_create(l_available(1:max_num_parcels), &
                                    win_size,                       &
                                    disp_unit,                      &
                                    MPI_INFO_NULL,                  &
                                    comm%world,                     &
                                    win_avail,                      &
                                    comm%err)
                call MPI_Win_create(l_merged(1:max_num_parcels),    &
                                    win_size,                       &
                                    disp_unit,                      &
                                    MPI_INFO_NULL,                  &
                                    comm%world,                     &
                                    win_merged,                     &
                                    comm%err)
            endif
        end subroutine nearest_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine nearest_deallocate
            if (allocated(nppc)) then
                call MPI_Win_free(win_leaf, comm%err)
                call mpi_check_for_error(&
                    "in MPI_Win_free of parcel_nearest::nearest_deallocate.")

                call MPI_Win_free(win_avail, comm%err)
                call mpi_check_for_error(&
                    "in MPI_Win_free of parcel_nearest::nearest_deallocate.")

                call MPI_Win_free(win_merged, comm%err)
                call mpi_check_for_error(&
                    "in MPI_Win_free of parcel_nearest::nearest_deallocate.")

                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
                deallocate(l_leaf)
                deallocate(l_available)
                deallocate(l_merged)
            endif

            if (allocated(rsma)) then
                deallocate(rsma)
                deallocate(pidsma)
                deallocate(midsma)
            endif
        end subroutine nearest_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! @param[out] isma indices of small parcels
        ! @param[out] iclo indices of close parcels
        ! @param[out] inva indices of small + invalid (due to sending to other processes)
        !             parcels
        ! @param[out] n_local_small the array size of isma and iclo
        ! @param[out] n_invalid number of invalid parcels
        ! @post
        !   - isma and inva must be sorted in ascending order
        !   - isma, inva and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(isma, iclo, inva, n_local_small, n_invalid)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer, allocatable, intent(out) :: inva(:)
            integer,              intent(out) :: n_local_small
            integer,              intent(out) :: n_invalid
            integer                           :: n_global_small
            integer                           :: ijk, n, k, j
            integer, allocatable              :: rclo(:)    ! MPI rank of closest parcel
            double precision, allocatable     :: dclo(:)    ! distance to closest parcel

            call start_timer(merge_nearest_timer)

            call nearest_allocate

            ! We must store the parcel index and the merge index *m*
            ! of each small parcel. We do not need to allocate the
            ! invalid buffer, therefore the second argument is .false.
            call allocate_parcel_id_buffers(2, .false.)

            !---------------------------------------------------------------------
            ! Initialise search:
            n_invalid = 0
            n_local_small = 0
            n_global_small = 0
            n_neighbour_small = 0
            n_remote_small = 0
            nppc = 0 !nppc(ijk) will contain the number of parcels in grid cell ijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels

                call parcel_to_local_cell_index(n)

                if (parcels%volume(n) < vmin) then
                    n_local_small = n_local_small + 1

                    ! If a small parcel is in a boundary cell, a duplicate must
                    ! be sent to the neighbour rank. This call checks if the parcel
                    ! must be sent and fills the send buffers.
                    call locate_parcel_in_boundary_cell(n_local_small, n)
                endif
            enddo

            call MPI_Allreduce(n_local_small,   &
                               n_global_small,  &
                               1,               &
                               MPI_INTEGER,     &
                               MPI_SUM,         &
                               comm%world,      &
                               comm%err)

            call mpi_check_for_error(&
                    "in MPI_Allreduce of parcel_nearest::find_nearest.")

            if (n_global_small == 0) then
                call nearest_deallocate
                call deallocate_parcel_id_buffers
                call stop_timer(merge_nearest_timer)
                return
            endif

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute, parcel index and the merge index of small parcels
            ! in boundary region to neighbours. We only need to send the position
            ! as this is the only attribute needed to figure out with whom a parcel
            ! might merge.
            call send_small_parcel_bndry_info

            ! There are 4 cases:
            !   - n_local_small = 0 and n_remote_small = 0 --> This rank has no small parcels.
            !   - n_local_small > 0 and n_remote_small = 0 --> This rank has only local small parcels.
            !   - n_local_small = 0 and n_remote_small > 0 --> This rank has only remote small parcels.
            !   - n_local_small > 0 and n_remote_small > 0 --> This rank has local and remote small parcels.
            if (n_local_small + n_remote_small == 0) then
                ! Although *this* rank has no small parcels and is therefore not involved
                ! in any merging, it must nonetheless call the tree resolving routine as
                ! there is global synchronisation taking place.
                call resolve_tree(isma, iclo, rclo, n_local_small)

                return
            endif

            ! allocate arrays
            allocate(isma(0:n_local_small + n_remote_small))
            allocate(inva(0:n_local_small + n_remote_small))
            allocate(iclo(n_local_small + n_remote_small))
            allocate(rclo(n_local_small + n_remote_small))
            allocate(dclo(n_local_small + n_remote_small))

            isma = -1
            iclo = -1
            rclo = -1
            dclo = product(extent)

            ! Find arrays kc1(ijk) & kc2(ijk) which indicate the parcels in grid cell ijk
            ! through n = node(k), for k = kc1(ijk), kc2(ijk):
            kc1(1) = 1
            do ijk = 1, box%halo_ncell-1
                kc1(ijk+1) = kc1(ijk) + nppc(ijk)
            enddo

            kc2 = kc1 - 1
            j = 0
            do n = 1, n_parcels + n_remote_small
                ijk = loca(n)
                k = kc2(ijk) + 1
                node(k) = n
                kc2(ijk) = k

                if (parcels%volume(n) < vmin) then
                    j = j + 1
                    isma(j) = n
                endif
            enddo

            !---------------------------------------------------------------------
            ! Determine locally closest parcel:
            call find_closest_parcel_locally(n_local_small, isma, iclo, rclo, dclo)

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(n_local_small, iclo, rclo, dclo)

            call stop_timer(merge_nearest_timer)

            !---------------------------------------------------------------------
            ! Figure out the mergers:
            call resolve_tree(isma, iclo, rclo, n_local_small)

            !---------------------------------------------------------------------
            ! Mark all entries of isma, iclo and rclo above n_local_small
            ! as invalid.
            do n = n_local_small+1, size(iclo)
                isma(n) = -1
                iclo(n) = -1
                rclo(n) = -1
            enddo

            isma(1:n_local_small) = pack(isma(1:), isma(1:) /= -1)
            iclo(1:n_local_small) = pack(iclo, iclo /= -1)
            rclo(1:n_local_small) = pack(rclo, rclo /= -1)

            call deallocate_parcel_id_buffers

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: It cannot happen that the closest parcel is a small parcel
            !       on another MPI rank that is sent elsewhere.
            call gather_remote_parcels(n_local_small, n_invalid, rclo, iclo, isma, inva)

            !------------------------------------------------------------------
            ! Sanity check: Indices of close parcels must be smaller equal to the
            ! number of local parcels:
            if (maxval(iclo(1:n_local_small)) > n_parcels) then
                call mpi_exit_on_error(&
                    "in parcel_nearest::find_nearest: Close parcel index out of range.")
            endif
            call nearest_deallocate

        end subroutine find_nearest

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! If a parcel communication happens across a periodic boundary, we must
        ! shift the parcel position such that it ends up in the sub-domain.
        subroutine apply_periodic_shift(pos, dir)
            double precision, intent(inout) :: pos(3)
            integer,          intent(in)    :: dir
            logical                         :: l_x_edge, l_y_edge
            integer                         :: ix, iy

            ix = nint(dxi(1) * (pos(1) - lower(1)))
            iy = nint(dxi(2) * (pos(2) - lower(2)))

            ! put <= and >= for safety reasons; it should actually always be ==
            l_x_edge = (ix <= 0) .or. (ix >= nx)
            l_y_edge = (iy <= 0) .or. (iy >= ny)

            select case (dir)
                case (MPI_NORTH)
                    pos(2) = merge(pos(2) + extent(2), pos(2), l_y_edge)
                case (MPI_SOUTH)
                    pos(2) = merge(pos(2) - extent(2), pos(2), l_y_edge)
                case (MPI_WEST)
                    pos(1) = merge(pos(1) - extent(1), pos(1), l_x_edge)
                case (MPI_EAST)
                    pos(1) = merge(pos(1) + extent(1), pos(1), l_x_edge)
                case (MPI_NORTHWEST)
                    pos(1) = merge(pos(1) - extent(1), pos(1), l_x_edge)
                    pos(2) = merge(pos(2) + extent(2), pos(2), l_y_edge)
                case (MPI_NORTHEAST)
                    pos(1) = merge(pos(1) + extent(1), pos(1), l_x_edge)
                    pos(2) = merge(pos(2) + extent(2), pos(2), l_y_edge)
                case (MPI_SOUTHWEST)
                    pos(1) = merge(pos(1) - extent(1), pos(1), l_x_edge)
                    pos(2) = merge(pos(2) - extent(2), pos(2), l_y_edge)
                case (MPI_SOUTHEAST)
                    pos(1) = merge(pos(1) + extent(1), pos(1), l_x_edge)
                    pos(2) = merge(pos(2) - extent(2), pos(2), l_y_edge)
                case default
                    call mpi_exit_on_error(&
                        "in parcel_nearest::apply_periodic_shift: No valid direction.")
            end select
        end subroutine apply_periodic_shift

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !!@pre Assumes a parcel is in the local domain including halo cells
        !      (in x and y).
        subroutine parcel_to_local_cell_index(n)
            integer, intent(in)  :: n
            integer              :: ix, iy, iz, ijk

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
        subroutine find_closest_parcel_locally(n_local_small, isma, iclo, rclo, dclo)
            integer,          intent(in)    :: n_local_small
            integer,          intent(inout) :: isma(0:)
            integer,          intent(inout) :: iclo(:)
            integer,          intent(inout) :: rclo(:)
            double precision, intent(inout) :: dclo(:)
            integer                         :: ix0, iy0, iz0, ijk, n, ix, iy, iz, m, k, ic, is
            double precision                :: delx, dely, delz, dsq, dsqmin, x_small, y_small, z_small

            ! SB: do not use temporary (j) index here, so we will be able to parallelise.
            ! Rather, stop if no nearest parcel found  in surrounding grid boxes
            do m = 1, n_local_small + n_remote_small
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
                    do iy = max(0, iy0-1), min(box%size(2)+1, iy0)
                        do ix = max(0, ix0-1), min(box%size(1)+1, ix0)
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
                rclo(m) = comm%rank
                dclo(m) = dsqmin
                l_merged(is) = .false.
                l_merged(ic) = .false.
            enddo

            !---------------------------------------------------------------------
            ! Update isma, iclo and rclo with indices of remote parcels:
            do m = 1, n_local_small + n_remote_small
                is = isma(m)
                ic = iclo(m)

                if ((is > n_parcels) .and. (ic > n_parcels)) then
                    ! A remote small parcel points to another remote small parcel. The remotes do not necessarily
                    ! need to be the same. Also, it can be a dual-link. As the same distance is evaluated on
                    ! the other two MPI ranks, *this* MPI rank must set the distance between the parcels to the
                    ! maximum value as otherwise the function "find_closest_parcel_globally" may think the
                    ! parcels belong to *this* MPI rank due to round-offs in the distance calculation. The
                    ! function "find_closest_parcel_globally" always sets "rclo" to the MPI source.
                    dclo(m) = huge(0.0d0) ! huge(x) returns the maximum value of this type
                else if (ic > n_parcels) then
                    ! A local small parcel points to a remote small parcel.
                    ! The index *ic* is larger than the local number of parcels
                    ! we must update *iclo(m)* and *rclo(m)* with the index of the parcel stored
                    ! on the other MPI rank.
                    iclo(m) = pidsma(ic)
                    rclo(m) = rsma(ic)
                endif
            enddo

            !------------------------------------------------------------------
            ! Sanity check: Indices of small parcels must be smaller equal to the
            ! number of local parcels:
            do m = 1, n_local_small
                if (isma(m) > n_parcels) then
                    call mpi_exit_on_error(&
                        'in in parcel_nearest::find_closest_parcel_locally: Small parcel index out of range.')
                endif
            enddo
        end subroutine find_closest_parcel_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! For all small parcels in boundary cells we must now determine the
        ! closest parcel globally (that is we must compare with the neighbour rank).
        ! For this purpose we send the information of the small parcels that we received
        ! back to the original ranks which then figure out the closest parcel.
        ! Note: The information about the received small parcels is stored in the last
        !       n_remote_small entries of isma, iclo, rclo and dlco
        subroutine find_closest_parcel_globally(n_local_small, iclo, rclo, dclo)
            integer,          intent(in)            :: n_local_small
            integer,          intent(inout)         :: iclo(:)
            integer,          intent(inout)         :: rclo(:)
            double precision, intent(inout)         :: dclo(:)
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: recv_size, send_size, buf_sizes(8)
            integer                                 :: tag, source, recv_count, n, l, i, m, k, j
            integer, parameter                      :: n_entries = 3

            buf_sizes = n_neighbour_small * n_entries
            call allocate_mpi_buffers(buf_sizes)

            !------------------------------------------------------------------
            ! Communicate with neighbours:
            j = 0
            do n = 1, 8

                ! we send the distance, the remote parcel index and the
                ! remote merge index *m* --> n_entries = 3.
                send_size = n_neighbour_small(n) * n_entries

                call get_mpi_buffer(n, send_buf)

                do l = 1, n_neighbour_small(n)
                    i = 1 + (l-1) * n_entries
                    k = j + n_parcels + l


                    ! merge index on *this* rank
                    m = n_local_small + j + l

                    send_buf(i)   = dble(midsma(k)) ! merge index on remote rank
                    send_buf(i+1) = dclo(m)         ! distance to closest parcel
                    send_buf(i+2) = dble(iclo(m))   ! parcel index of closest parcel
                enddo

                j = j + n_neighbour_small(n)

                tag = small_recv_order(n)

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(tag)%rank,    &
                               NEIGHBOUR_TAG(tag),      &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error(&
                    "in MPI_Isend of parcel_nearest::find_closest_parcel_globally.")

            enddo

            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::find_closest_parcel_globally: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                call mpi_check_for_error(&
                    "in MPI_Recv of parcel_nearest::find_closest_parcel_globally.")

                do l = 1, recv_count
                    i = 1 + (l-1) * n_entries
                    m = nint(recv_buf(i))
                    if (dclo(m) > recv_buf(i+1)) then
                        ! the local closest distance is
                        ! larger; we must use the remote parcel and
                        ! therefore update the rclo entry.
                        dclo(m) = recv_buf(i+1)
                        iclo(m) = nint(recv_buf(i+2))
                        rclo(m) = source
                    endif
                enddo

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            call mpi_check_for_error(&
                "in MPI_Waitall of parcel_nearest::find_closest_parcel_globally.")

            call deallocate_mpi_buffers

        end subroutine find_closest_parcel_globally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! https://github.com/mpi-forum/mpi-forum-historic/issues/413
        ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node294.htm
        ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node279.htm
        subroutine resolve_tree(isma, iclo, rclo, n_local_small)
            integer, intent(inout)         :: isma(0:)
            integer, intent(inout)         :: iclo(:)
            integer, intent(inout)         :: rclo(:)
            integer, intent(inout)         :: n_local_small
            integer                        :: ic, rc, is, m, j
            logical                        :: l_helper
            integer(KIND=MPI_ADDRESS_KIND) :: offset
            logical                        :: l_continue_iteration, l_do_merge(n_local_small)
            logical                        :: l_isolated_dual_link(n_local_small)

            call start_timer(merge_tree_resolve_timer)

            ! First, iterative, stage
            l_continue_iteration = .true.

            do while (l_continue_iteration)
                l_continue_iteration = .false.
                ! reset relevant properties for candidate mergers

                ! synchronize the private and public window copies
                call MPI_Win_sync(win_merged, comm%err)
                call MPI_Win_sync(win_avail, comm%err)
                call MPI_Win_sync(win_leaf, comm%err)

                do m = 1, n_local_small
                    is = isma(m)
                    ! only consider links that still may be merging
                    ! reset relevant properties
                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)
                        l_leaf(is) = .true.

                        if (rc == comm%rank) then
                            l_available(ic) = .true.
                        else
                            call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)
                            !     MPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                            !         target_disp, target_count, target_datatype, win, ierror)
                            !     TYPE(*), DIMENSION(..), INTENT(IN), ASYNCHRONOUS :: origin_addr
                            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
                            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
                            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
                            !     TYPE(MPI_Win), INTENT(IN) :: win
                            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                            l_helper = .true.
                            offset = ic - 1 ! starts at 0
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         offset,            &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)
                            call mpi_check_for_error(&
                                "in MPI_Put of parcel_nearest::resolve_tree.")

                            call MPI_Win_unlock(rc, win_avail, comm%err)
                        endif
                    endif
                enddo

                call MPI_Barrier(comm%world, comm%err)

                ! synchronize the private and public window copies
                call MPI_Win_sync(win_avail, comm%err)

                ! determine leaf parcels
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        if (rc == comm%rank) then
                            l_leaf(ic) = .false.
                        else
                            call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_leaf, comm%err)
                            l_helper = .false.
                            offset = ic - 1
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         offset,            &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         win_leaf,          &
                                         comm%err)
                            call mpi_check_for_error(&
                                "in MPI_Put of parcel_nearest::resolve_tree.")
                            call MPI_Win_unlock(rc, win_leaf, comm%err)
                        endif
                    endif
                enddo

                call MPI_Barrier(comm%world, comm%err)

                ! synchronize the private and public window copies
                call MPI_Win_sync(win_leaf, comm%err)

                ! filter out parcels that are "unavailable" for merging
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. l_merged(is)) then
                        if (.not. l_leaf(is)) then
                            ic = iclo(m)
                            rc = rclo(m)

                            if (rc == comm%rank) then
                                l_available(ic) = .false.
                            else
                                call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)
                                l_helper = .false.
                                offset = ic - 1
                                call MPI_Put(l_helper,          &
                                             1,                 &
                                             MPI_LOGICAL,       &
                                             rc,                &
                                             offset,            &
                                             1,                 &
                                             MPI_LOGICAL,       &
                                             win_avail,         &
                                             comm%err)
                                call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")
                                call MPI_Win_unlock(rc, win_avail, comm%err)
                            endif
                        endif
                    endif
                enddo

                call MPI_Barrier(comm%world, comm%err)

                ! synchronize the private and public window copies
                call MPI_Win_sync(win_avail, comm%err)

                ! identify mergers in this iteration
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        l_helper = .false.
                        if (rc == comm%rank) then
                            l_helper = l_available(ic)
                        else
                            call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)
                            !     MPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                            !         target_disp, target_count, target_datatype, win, ierror)
                            !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: origin_addr
                            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
                            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
                            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
                            !     TYPE(MPI_Win), INTENT(IN) :: win
                            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                            offset = ic - 1
                            call MPI_Get(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         offset,            &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)
                            call mpi_check_for_error(&
                                    "in MPI_Get of parcel_nearest::resolve_tree.")

                            call MPI_Win_unlock(rc, win_avail, comm%err)
                        endif

                        if (l_leaf(is) .and. l_helper) then
                            l_continue_iteration = .true. ! merger means continue iteration
                            l_merged(is) = .true.

                            if (rc == comm%rank) then
                                l_merged(ic) = .true.
                            else
                                call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_merged, comm%err)

                                l_helper = .true.
                                offset = ic - 1
                                call MPI_Put(l_helper,          &
                                            1,                  &
                                            MPI_LOGICAL,        &
                                            rc,                 &
                                            offset,             &
                                            1,                  &
                                            MPI_LOGICAL,        &
                                            win_merged,         &
                                            comm%err)
                                call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")
                                call MPI_Win_unlock(rc, win_merged, comm%err)
                            endif
                        endif
                    endif
                enddo

                ! Performance improvement: We actually only need to synchronize with neighbours
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   l_continue_iteration,    &
                                   1,                       &
                                   MPI_LOGICAL,             &
                                   MPI_LOR,                 &
                                   comm%world,              &
                                   comm%err)

                call mpi_check_for_error(&
                    "in MPI_Allreduce of parcel_nearest::resolve_tree.")
            enddo

            call MPI_Barrier(comm%world, comm%err)

            ! synchronize the private and public window copies
            call MPI_Win_sync(win_merged, comm%err)
            call MPI_Win_sync(win_avail, comm%err)
            call MPI_Win_sync(win_leaf, comm%err)


            ! Second stage, related to dual links
            do m = 1, n_local_small
                is = isma(m)

                if (.not. l_merged(is)) then
                    if (l_leaf(is)) then ! set in last iteration of stage 1
                        ic = iclo(m)
                        rc = rclo(m)

                        if (rc == comm%rank) then
                            l_available(ic) = .true.
                        else
                            call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)
                            l_helper = .true.
                            offset = ic - 1
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         offset,            &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)

                            call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")

                            call MPI_Win_unlock(rc, win_avail, comm%err)
                        endif
                    endif
                endif
            enddo

            call MPI_Barrier(comm%world, comm%err)

            ! synchronize the private and public window copies
            call MPI_Win_sync(win_avail, comm%err)

            ! Second stage
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                l_do_merge(m) = .false.
                l_isolated_dual_link(m) = .false.

                if (l_merged(is) .and. l_leaf(is)) then
                    ! previously identified mergers: keep
                    l_do_merge(m) = .true.
                    !----------------------------------------------------------
                    ! begin of sanity check
                    ! After first stage mergers parcel cannot be both initiator
                    ! and receiver in stage 1
                    l_helper = .false.
                    if (rc == comm%rank) then
                        l_helper = l_leaf(ic)
                    else
                        call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_leaf, comm%err)
                        offset = ic - 1
                        call MPI_Get(l_helper,          &
                                     1,                 &
                                     MPI_LOGICAL,       &
                                     rc,                &
                                     offset,            &
                                     1,                 &
                                     MPI_LOGICAL,       &
                                     win_leaf,          &
                                     comm%err)
                        call mpi_check_for_error(&
                            "in MPI_Get of parcel_nearest::resolve_tree.")

                        call MPI_Win_unlock(rc, win_leaf, comm%err)
                    endif
                    if (l_helper) then
                        call mpi_exit_on_error(&
                            'in parcel_nearest::resolve_tree: First stage error')
                    endif
                    ! end of sanity check
                    !----------------------------------------------------------

                elseif (.not. l_merged(is)) then
                    if (l_leaf(is)) then
                        ! links from leafs
                        l_do_merge(m) = .true.
                    elseif (.not. l_available(is)) then
                        ! Above means parcels that have been made 'available' do not keep outgoing links

                        if (rc == comm%rank) then
                            l_helper = l_available(ic)
                        else
                            call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)

                            offset = ic - 1
                            call MPI_Get(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         offset,            &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)

                            call mpi_check_for_error(&
                                "in MPI_Get of parcel_nearest::resolve_tree.")

                            call MPI_Win_unlock(rc, win_avail, comm%err)
                        endif

                        if (l_helper) then
                            ! merge this parcel into ic along with the leaf parcels
                            l_do_merge(m) = .true.

                        else
                            l_isolated_dual_link(m) = .true.
                            ! isolated dual link
                            ! Don't keep current link
                            ! But make small parcel available so other parcel can merge with it
                            ! THIS NEEDS THINKING ABOUT A PARALLEL IMPLEMENTATION
                            ! This could be based on the other parcel being outside the domain
                            ! And a "processor order"
                            if (comm%rank <= rc) then
                                ! The MPI rank with lower number makes its parcel
                                ! available.
                                l_available(is) = .true.
                            endif
                        endif
                    endif
                endif
            enddo

            call MPI_Barrier(comm%world, comm%err)

            ! synchronize the private and public window copies
            call MPI_Win_sync(win_avail, comm%err)


            !------------------------------------------------------
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)

                if ((l_do_merge(m) .eqv. .false.) .and. l_isolated_dual_link(m)) then
                    ! isolated dual link

                    if (rc == comm%rank) then
                        l_helper = l_available(ic)
                    else
                        call MPI_Win_lock(MPI_LOCK_SHARED, rc, 0, win_avail, comm%err)
                        offset = ic - 1
                        call MPI_Get(l_helper,          &
                                     1,                 &
                                     MPI_LOGICAL,       &
                                     rc,                &
                                     offset,            &
                                     1,                 &
                                     MPI_LOGICAL,       &
                                     win_avail,         &
                                     comm%err)

                        call mpi_check_for_error(&
                            "in MPI_Get of parcel_nearest::resolve_tree.")

                        call MPI_Win_unlock(rc, win_avail, comm%err)
                    endif

                    if (l_helper) then
                        ! merge this parcel into ic along with the leaf parcels
                        l_do_merge(m) = .true.
!                     else
!                       ! Dual link is resolved on other rank
                    endif
                endif
                !------------------------------------------------------
            enddo

            call MPI_Barrier(comm%world, comm%err)

            ! synchronize the private and public window copies
            call MPI_Win_sync(win_avail, comm%err)

            j = 0
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                if (l_do_merge(m)) then
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic
                    rclo(j) = rc
                endif
            enddo
            n_local_small = j

            call stop_timer(merge_tree_resolve_timer)
        end subroutine resolve_tree

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This routine fills the *pid* buffers with 2 entries per small parcel.
        ! The first entry is the local parcel index in the parcel container and
        ! the second entry is the merge index (usually accessed with *m*).
        subroutine locate_parcel_in_boundary_cell(m, n)
            integer, intent(in) :: m, n
            integer             :: k, ix, iy

            ! nearest global grid point
            ix = nint(dxi(1) * (parcels%position(1, n) - lower(1)))
            iy = nint(dxi(2) * (parcels%position(2, n) - lower(2)))

            ! check lower x-direction
            if (ix == box%lo(1)) then

                if (iy == box%lo(2)) then
                    ! parcel in southwest corner with
                    ! neighbours: west, south and southwest
                    k = n_parcel_sends(MPI_SOUTHWEST) + 1
                    southwest_pid(2*k-1) = n
                    southwest_pid(2*k) = m
                    n_parcel_sends(MPI_SOUTHWEST) = k

                else if (iy == box%hi(2)+1) then
                    ! parcel in northwest corner with
                    ! neighbours: west, north and northwest
                    k = n_parcel_sends(MPI_NORTHWEST) + 1
                    northwest_pid(2*k-1) = n
                    northwest_pid(2*k) = m
                    n_parcel_sends(MPI_NORTHWEST) = k
                endif

                ! neighbour: west
                k = n_parcel_sends(MPI_WEST) + 1
                west_pid(2*k-1) = n
                west_pid(2*k) = m
                n_parcel_sends(MPI_WEST) = k
            endif

            ! check upper x-direction (use >= as a parcel might be directly on the cell edge)
            if (ix == box%hi(1)+1) then

                if (iy == box%lo(2)) then
                    ! parcel in southeast corner with
                    ! neighbours: east, south and southeast
                    k = n_parcel_sends(MPI_SOUTHEAST) + 1
                    southeast_pid(2*k-1) = n
                    southeast_pid(2*k) = m
                    n_parcel_sends(MPI_SOUTHEAST) = k

                else if (iy == box%hi(2)+1) then
                    ! parcel in northeast corner with
                    ! neighbours: east, north and northeast
                    k = n_parcel_sends(MPI_NORTHEAST) + 1
                    northeast_pid(2*k-1) = n
                    northeast_pid(2*k) = m
                    n_parcel_sends(MPI_NORTHEAST) = k
                endif

                ! neighbours: east
                k = n_parcel_sends(MPI_EAST) + 1
                east_pid(2*k-1) = n
                east_pid(2*k) = m
                n_parcel_sends(MPI_EAST) = k
            endif

            if (iy == box%lo(2)) then
                k = n_parcel_sends(MPI_SOUTH) + 1
                south_pid(2*k-1) = n
                south_pid(2*k) = m
                n_parcel_sends(MPI_SOUTH) = k
            endif

            if (iy == box%hi(2)+1) then
                k = n_parcel_sends(MPI_NORTH) + 1
                north_pid(2*k-1) = n
                north_pid(2*k) = m
                n_parcel_sends(MPI_NORTH) = k
            endif
        end subroutine locate_parcel_in_boundary_cell

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Send position of small parcels in boundary region, their local index
        ! in the parcel container and their merge index (*m*) to neighbours
        subroutine send_small_parcel_bndry_info
            integer,          dimension(:), pointer :: send_ptr
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: recv_size, send_size
            integer                                 :: tag, source, recv_count, n, i ,j, l, m, pid, k
            integer, parameter                      :: n_entries = 5
            integer, allocatable                    :: rtmp(:), pidtmp(:), midtmp(:)
            ! rtmp: MPI rank remote small parcel belongs to
            ! pidtmp: parcel index of remote parcel (on the owning rank)
            ! midtmp: m index of remote parcel (on the owning rank)

            n_neighbour_small = 0

            ! dummy allocate; will be updated on the fly
            allocate(pidsma(0))
            allocate(rsma(0))
            allocate(midsma(0))

            ! we only send parcel position, parcel index and its m index (of isma(m))
            call allocate_parcel_buffers(n_entries)

            do n = 1, 8

                ! Each entry of the buffer send_ptr points to has 2 entries
                ! per small parcel: 1. local parcel index; 2. merge index *m*
                call get_parcel_buffer_ptr(n, send_ptr, send_buf)

                send_size = n_parcel_sends(n) * n_entries

                if (n_parcel_sends(n) > 0) then
                    ! pack parcel position and parcel index to send buffer
                    do l = 1, n_parcel_sends(n)
                        i = 1 + (l-1) * n_entries
                        pid = send_ptr(2*l-1)
                        m   = send_ptr(2*l)
                        send_buf(i:i+2) = parcels%position(:, pid)
                        send_buf(i+3) = dble(pid)
                        send_buf(i+4) = dble(m)
                    enddo
                endif

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error(&
                    "in MPI_Isend of parcel_nearest::send_small_parcel_bndry_info.")
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                call mpi_check_for_error(&
                    "in MPI_Recv of parcel_nearest::send_small_parcel_bndry_info.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::send_small_parcel_bndry_info: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                i = n_parcels+1
                j = n_parcels+1+size(rsma) + recv_count
                allocate(rtmp(i:j))
                allocate(pidtmp(i:j))
                allocate(midtmp(i:j))

                ! copy old over
                if (size(rsma) > 0) then
                    rtmp(n_parcels+1:n_parcels+size(rsma)) = rsma
                    pidtmp(n_parcels+1:n_parcels+size(pidsma)) = pidsma
                    midtmp(n_parcels+1:n_parcels+size(midsma)) = midsma
                endif

                deallocate(pidsma)
                deallocate(rsma)
                deallocate(midsma)

                if (recv_count > 0) then
                    ! unpack parcel position and parcel index to recv buffer
                    do l = 1, recv_count
                        i = 1 + (l-1) * n_entries
                        k = sum(n_neighbour_small) + n_parcels + l
                        parcels%position(:, k) = recv_buf(i:i+2)
                        parcels%volume(k) = zero    ! set to zero as each parcel is small
                        pidtmp(k) = nint(recv_buf(i+3))
                        rtmp(k) = source
                        midtmp(k) = nint(recv_buf(i+4))

                        !------------------------------------------------------
                        ! We must also assign incoming small parcels to cells
                        ! Note: We must apply a shift for parcels communicated
                        !       across a periodic boundary.
                        call apply_periodic_shift(parcels%position(:, k), tag)

                        ! Now we can safely assign the local cell index:
                        call parcel_to_local_cell_index(k)
                    enddo
                    n_neighbour_small(n) = n_neighbour_small(n) + recv_count
                endif

                small_recv_order(n) = tag

                ! *move_alloc* deallocates pidtmp, rtmp and midtmp
                call move_alloc(pidtmp, pidsma)
                call move_alloc(rtmp, rsma)
                call move_alloc(midtmp, midsma)

                deallocate(recv_buf)
            enddo

            n_remote_small = sum(n_neighbour_small)

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            call mpi_check_for_error(&
                "in MPI_Waitall of parcel_nearest::send_small_parcel_bndry_info.")


            call deallocate_parcel_buffers

        end subroutine send_small_parcel_bndry_info

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine gather_remote_parcels(n_local_small, n_invalid, rclo, iclo, isma, inva)
            integer,                     intent(inout) :: n_local_small
            integer,                     intent(inout) :: n_invalid
            integer,                     intent(inout) :: rclo(:)       ! MPI rank of closest parcel
            integer,                     intent(inout) :: iclo(:)
            integer,                     intent(inout) :: isma(0:)
            integer,                     intent(inout) :: inva(0:)
            double precision, dimension(:), pointer    :: send_buf
            double precision, allocatable              :: recv_buf(:)
            type(MPI_Request)                          :: requests(8)
            type(MPI_Status)                           :: recv_status, send_statuses(8)
            integer                                    :: tag, source, recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(n_par_attrib)
            integer                                    :: m, rc, ic, is, n, i, j, k, iv
            integer                                    :: n_entries
            integer                                    :: n_bytes
            integer, asynchronous                      :: n_parcel_recvs(8)
            type(MPI_Win)                              :: win_neighbour
            integer(KIND=MPI_ADDRESS_KIND)             :: win_size, offset
            integer                                    :: n_registered(8)

            !------------------------------------------------------------------
            ! Figure out how many parcels we send:
            n_parcel_sends = 0

            do m = 1, n_local_small
                rc = rclo(m)
                if (comm%rank /= rc) then
                    n = get_neighbour_from_rank(rc)
                    n_parcel_sends(n) = n_parcel_sends(n) + 1
                endif
            enddo

            n_invalid = sum(n_parcel_sends)

            !------------------------------------------------------------------
            ! We must send all parcel attributes (n_par_attrib) plus
            ! the index of the close parcel ic (1)
            n_entries = n_par_attrib + 1
            n_registered = n_parcel_sends * n_entries
            call allocate_mpi_buffers(n_registered)

            !------------------------------------------------------------------
            ! Figure out how many parcels we receive from our neighbours
            ! and fill all buffers:
            call MPI_Sizeof(m, n_bytes, comm%err)
            win_size = 8 * n_bytes
            call MPI_Win_create(n_parcel_recvs,   &
                                win_size,         &
                                n_bytes,          &
                                MPI_INFO_NULL,    &
                                comm%world,       &
                                win_neighbour,    &
                                comm%err)

            n_parcel_recvs = 0
            n_registered = 0
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                inva(m) = is ! copy initial isma to inva for parcel deletion
                if (comm%rank /= rc) then
                    ! The closest parcel to this small parcel *is*
                    ! is on another MPI rank. We must send this parcel
                    ! to that rank.

                    n = get_neighbour_from_rank(rc)

                    call get_mpi_buffer(n, send_buf)

                    n_registered(n) = n_registered(n) + 1

                    call parcel_serialize(is, buffer)
                    j = n_entries * n_registered(n)
                    i = j - n_entries + 1
                    send_buf(i:j-1) = buffer
                    send_buf(j) = dble(ic)

                    ! Mark this as invalid
                    isma(m) = -1
                    iclo(m) = -1
                    rclo(m) = -1
                endif
            enddo

            ! We must now remove all invalid entries in isma and
            ! iclo and also update the value of n_local_small:
            iv = n_local_small  ! we need to keep the original value for filling *inva*
            n_local_small = count(isma(1:) /= -1)
            isma(1:n_local_small) = pack(isma(1:), isma(1:) /= -1)
            iclo(1:n_local_small) = pack(iclo, iclo /= -1)
            rclo(1:n_local_small) = pack(rclo, rclo /= -1)

            call MPI_Win_lock_all(0, win_neighbour, comm%err)

            do n = 1, 8
                offset = n - 1
                call MPI_Put(n_parcel_sends(n),    &
                             1,                    &
                             MPI_INTEGER,          &
                             neighbours(n)%rank,   &
                             offset,               &
                             1,                    &
                             MPI_INTEGER,          &
                             win_neighbour,        &
                             comm%err)
            enddo

            call mpi_check_for_error(&
                "in MPI_Put of parcel_nearest::resolve_tree.")

            call MPI_Win_unlock_all(win_neighbour, comm%err)

            call MPI_Win_sync(win_neighbour, comm%err)

            call MPI_Win_free(win_neighbour, comm%err)

            !------------------------------------------------------------------
            ! Communicate parcels:
            do n = 1, 8
                call get_mpi_buffer(n, send_buf)

                send_size = n_parcel_sends(n) * n_entries

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error(&
                    "in MPI_Isend of parcel_nearest::gather_remote_parcels.")
            enddo

            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)

                call mpi_check_for_error(&
                    "in MPI_Recv of parcel_nearest::gather_remote_parcels.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::gather_remote_parcels: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                ! Set the current value of the *m* index to
                ! its current maximum value (excluding small
                ! parcel received earlier from other MPI ranks.
                m = n_local_small
                do k = 1, recv_count
                    ! We receive a small parcel;
                    ! append it to the container with index
                    ! "n_parcels+1"
                    n_parcels = n_parcels + 1
                    j = n_entries * k
                    i = j - n_entries + 1
                    buffer = recv_buf(i:j-1)
                    call parcel_deserialize(n_parcels, buffer)

                    ! Add the small parcel to isma and inva, and
                    ! its closest parcel to iclo
                    m = m + 1
                    iclo(m) = nint(recv_buf(j))
                    isma(m) = n_parcels
                    iv = iv + 1
                    inva(iv) = n_parcels
                enddo
                ! The last value of *m* is our new number of
                ! small parcels. However, this still includes
                ! invalid entries in isma and iclo that we will
                ! remove next.
                n_local_small = m

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            call mpi_check_for_error(&
                "in MPI_Waitall of parcel_nearest::gather_remote_parcels.")

            call deallocate_mpi_buffers

#ifndef NDEBUG
            n = n_parcels - n_invalid
            call mpi_blocking_reduce(n, MPI_SUM)
            if ((comm%rank == comm%master) .and. (.not. n == n_total_parcels)) then
                call mpi_exit_on_error(&
                    "in parcel_nearest::gather_remote_parcels: We lost parcels.")
            endif
#endif
        end subroutine gather_remote_parcels

end module parcel_nearest
