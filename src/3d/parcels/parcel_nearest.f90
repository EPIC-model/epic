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
#ifdef ENABLE_VERBOSE
    use options, only : verbose, output
#endif
    use constants, only : zero
    use parcels_mod
    use parcel_ops, only : get_delx &
                         , get_dely
    use fields, only : get_index
    use parameters, only : dx, dxi, vcell, hli, lower, extent       &
                         , ncell, nx, ny, nz, max_num_parcels
    use options, only : parcel
    use mpi_timer, only : start_timer, stop_timer, timings
    use mpi_environment
    use mpi_layout
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error          &
                        , mpi_check_for_error        &
                        , mpi_check_for_message      &
                        , mpi_stop                   &
                        , mpi_check_rma_window_model
    use iso_c_binding, only : c_ptr, c_f_pointer
    use parcel_nearest_mpi, only : tree_t
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
#ifndef NDEBUG
    use datatypes, only : int64
#endif
    implicit none

    integer:: merge_nearest_timer, merge_tree_resolve_timer
    integer:: nearest_allreduce_timer
    integer:: nearest_exchange_timer

    private

    type nearest_type
        integer, allocatable :: nppc(:), kc1(:), kc2(:)
        integer, allocatable :: loca(:)
        integer, allocatable :: node(:)

        integer, allocatable :: rsma(:)               ! rank of small parcel received (accessed with *n*)
        integer, allocatable :: pidsma(:)             ! index of remote small parcel (accessed with *n*)
        integer, allocatable :: midsma(:)             ! *m* in *isma(m)*

        contains
            procedure :: alloc => nearest_allocate
            procedure :: dealloc => nearest_deallocate
            procedure :: parcel_to_local_cell_index => nearest_parcel_to_local_cell_index
            procedure :: set_lbound => nearest_set_lbound
            procedure :: set_ubound => nearest_set_ubound

    end type nearest_type

    type(nearest_type) :: near

    type(tree_t) :: tree

    integer              :: n_neighbour_small(8)  ! number of small parcels received

    type(communicator) :: subcomm

#ifdef ENABLE_VERBOSE
    double precision :: simtime
#endif

#ifndef NDEBUG
    ! Small remote parcels must be sent back to the original
    ! MPI ranks in the same order as they were received. These
    ! two arrays verify the correctness.
    integer :: small_recv_order(8)
    integer :: small_recv_count(8)
#endif

    public :: find_nearest                      &
            , merge_nearest_timer               &
            , merge_tree_resolve_timer          &
            , update_remote_indices             &
            , locate_parcel_in_boundary_cell    &
            , send_small_parcel_bndry_info      &
            , find_closest_parcel_globally      &
            , nearest_allreduce_timer           &
#ifdef ENABLE_VERBOSE
            , simtime                           &
#endif
            , nearest_exchange_timer            &
            , near

    contains

        subroutine nearest_allocate(this, max_num_parcels)
            class(nearest_type), intent(inout) :: this
            integer,             intent(in)    :: max_num_parcels

            if (.not. allocated(this%nppc)) then
                allocate(this%nppc(box%halo_ncell))
                allocate(this%kc1(box%halo_ncell))
                allocate(this%kc2(box%halo_ncell))
                allocate(this%loca(max_num_parcels))
                allocate(this%node(max_num_parcels))
            endif
        end subroutine nearest_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine nearest_deallocate(this)
            class(nearest_type), intent(inout) :: this

            if (allocated(this%nppc)) then
                deallocate(this%nppc)
                deallocate(this%kc1)
                deallocate(this%kc2)
                deallocate(this%loca)
                deallocate(this%node)
            endif

            if (allocated(this%rsma)) then
                deallocate(this%rsma)
                deallocate(this%pidsma)
                deallocate(this%midsma)
            endif

        end subroutine nearest_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! If parcels are on the periodic edge, they might end up on the wrong side.
        ! We must correct this by shifting the parcel position such that it ends
        ! up in the sub-domain.
        subroutine handle_periodic_edge_parcels(pos)
            double precision, intent(inout) :: pos(:)
            integer                         :: ix, iy, iz

            ! global grid index in x and y
            call get_index(pos, ix, iy, iz)

            ! apply a periodic shift if not inside sub-domain
            ! note: we subtract 1 from box%hhi(1:2) as we have 2 halo grid points
            ! on the upper side
            pos(1) = merge(pos(1) + extent(1), pos(1), ix < box%hlo(1))
            pos(1) = merge(pos(1) - extent(1), pos(1), ix > box%hhi(1) - 1)

            pos(2) = merge(pos(2) + extent(2), pos(2), iy < box%hlo(2))
            pos(2) = merge(pos(2) - extent(2), pos(2), iy > box%hhi(2) - 1)

        end subroutine handle_periodic_edge_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !@pre Assumes a parcel is in the local domain including halo cells
        !      (in x and y).
        subroutine nearest_parcel_to_local_cell_index(this, pcont, n)
            class(nearest_type), intent(inout) :: this
            class(pc_type),      intent(inout) :: pcont
            integer,             intent(in)    :: n
            integer                            :: ix, iy, iz, ijk

            call handle_periodic_edge_parcels(pcont%position(:, n))

            ix =     int(dxi(1) * (pcont%position(1, n) - box%halo_lower(1)))
            iy =     int(dxi(2) * (pcont%position(2, n) - box%halo_lower(2)))
            iz = min(int(dxi(3) * (pcont%position(3, n) - box%lower(3))), nz-1)

            ! Cell index of parcel:
            !   This runs from 1 to halo_ncell where
            !   halo_ncell includes halo cells
            ijk = 1 + ix + box%halo_size(1) * iy + box%halo_size(1) * box%halo_size(2) * iz

#ifndef NEDBUG
            if (ijk < 1) then
               call mpi_exit_on_error(&
                        'in nearest_parcel_to_local_cell_index: Parcel not in local domain.')
            endif
#endif

            ! Accumulate number of parcels in this grid cell:
            this%nppc(ijk) = this%nppc(ijk) + 1

            ! Store grid cell that this parcel is in:
            this%loca(n) = ijk

        end subroutine nearest_parcel_to_local_cell_index

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Find arrays kc1(ijk) & kc2(ijk) which indicate the parcels in grid cell ijk
        ! through n = node(k), for k = kc1(ijk), kc2(ijk):
        subroutine nearest_set_lbound(this)
            class(nearest_type), intent(inout) :: this
            integer                            :: ijk

            this%kc1(1) = 1
            do ijk = 1, box%halo_ncell-1
                this%kc1(ijk+1) = this%kc1(ijk) + this%nppc(ijk)
            enddo

            this%kc2 = this%kc1 - 1

        end subroutine nearest_set_lbound

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine nearest_set_ubound(this, n)
            class(nearest_type), intent(inout) :: this
            integer,             intent(in)    :: n
            integer                            :: ijk, k

            ijk = this%loca(n)
            k = this%kc2(ijk) + 1
            this%node(k) = n
            this%kc2(ijk) = k

        end subroutine nearest_set_ubound

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
        subroutine find_nearest(pcont, isma, iclo, inva, n_local_small, n_invalid)
            class(pc_type),       intent(inout) :: pcont
            integer, allocatable, intent(out)   :: isma(:)
            integer, allocatable, intent(out)   :: iclo(:)
            integer, allocatable, intent(out)   :: inva(:)
            integer,              intent(out)   :: n_local_small
            integer,              intent(out)   :: n_invalid
            integer                             :: n_global_small, color
            integer                             :: n_remote_small        ! sum(n_neighbour_small)
            integer                             :: n, j
            integer, allocatable                :: rclo(:)    ! MPI rank of closest parcel
            double precision, allocatable       :: dclo(:)    ! distance to closest parcel
            logical                             :: l_no_small ! if *this* rank has no local and no remote small
                                                              ! parcels
#ifdef ENABLE_VERBOSE
            logical                           :: l_exist
            character(512)                    :: fname
#endif

            call start_timer(merge_nearest_timer)

            call near%alloc(pcont%max_num)

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

            call tree%initialise(max_num_parcels)

            near%nppc = 0 !nppc(ijk) will contain the number of parcels in grid cell ijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, pcont%local_num

                call near%parcel_to_local_cell_index(pcont, n)

                if (pcont%is_small(n)) then
                    n_local_small = n_local_small + 1

                    ! If a small parcel is in a boundary cell, a duplicate must
                    ! be sent to the neighbour rank. This call checks if the parcel
                    ! must be sent and fills the send buffers.
                    call locate_parcel_in_boundary_cell(pcont, n_local_small, n)
                endif
            enddo


            !------------------------------------------------------------------
            call MPI_Allreduce(n_local_small,   &
                               n_global_small,  &
                               1,               &
                               MPI_INTEGER,     &
                               MPI_SUM,         &
                               world%comm,      &
                               world%err)

            call mpi_check_for_error(world, &
                    "in MPI_Allreduce of parcel_nearest::find_nearest.")

            if (n_global_small == 0) then
                call near%dealloc
                call deallocate_parcel_id_buffers
                call tree%finalise
                call stop_timer(merge_nearest_timer)
                return
            endif

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute, parcel index and the merge index of small parcels
            ! in boundary region to neighbours. We only need to send the position
            ! as this is the only attribute needed to figure out with whom a parcel
            ! might merge.
            call send_small_parcel_bndry_info(pcont, n_remote_small)

            ! There are 4 cases:
            !   - n_local_small = 0 and n_remote_small = 0 --> This rank has no small parcels.
            !   - n_local_small > 0 and n_remote_small = 0 --> This rank has only local small parcels.
            !   - n_local_small = 0 and n_remote_small > 0 --> This rank has only remote small parcels.
            !   - n_local_small > 0 and n_remote_small > 0 --> This rank has local and remote small parcels.

            ! Although *this* rank has no small parcels and is therefore not involved
            ! in any merging, it must nonetheless call the tree resolving and parcel gathering
            ! routine as there is global synchronisation taking place.
            l_no_small = (n_local_small + n_remote_small == 0)

            !------------------------------------------------------------------
            ! Create subcommunicator:
            ! Only MPI ranks that have small parcels or received remote small parcels
            ! must be part of the communicator.

            ! Ensure the communicator is freed first.
            if (subcomm%comm /= MPI_COMM_NULL) then
                call MPI_Comm_free(subcomm%comm, subcomm%err)
                call mpi_check_for_error(subcomm, &
                        "in MPI_Comm_free of parcel_nearest::find_nearest.")
            endif

            ! Each MPI process must know if it is part of the subcommunicator or not.
            ! All MPI ranks that have small parcels or received small parcels from neighbouring
            ! MPI ranks must be part of the communicator.
            color = MPI_UNDEFINED
            if (.not. l_no_small) then
                color = 0  ! any non-negative number is fine
            endif

            call MPI_Comm_split(comm=cart%comm,         &
                                color=color,            &
                                key=cart%rank,          &  ! key controls the ordering of the processes
                                newcomm=subcomm%comm,   &
                                ierror=cart%err)

            if (subcomm%comm /= MPI_COMM_NULL) then
                ! The following two calls are not necessary, but we do for good practice.
                call MPI_Comm_size(subcomm%comm, subcomm%size, subcomm%err)
                call MPI_Comm_rank(subcomm%comm, subcomm%rank, subcomm%err)
                subcomm%root = 0

#ifdef ENABLE_VERBOSE
                if (verbose .and. (subcomm%rank == subcomm%root)) then
                    fname = trim(output%basename) // '_nearest_subcomm.asc'
                    inquire(file=trim(fname), exist=l_exist)
                    if (l_exist) then
                        open(unit=1236, file=trim(fname), status='old', position='append')
                    else
                        open(unit=1236, file=trim(fname), status='replace')
                        write(1236, *) '  # time            subcomm%size    percentage (%)'
                    endif
                    write(1236, *) simtime, subcomm%size, subcomm%size / dble(world%size) * 100.d0
                    close(1236)
                endif
#endif
            endif

            if (.not. l_no_small) then
                ! allocate arrays
                allocate(isma(0:n_local_small+n_remote_small))
                allocate(inva(0:n_local_small+n_remote_small))
                allocate(iclo(n_local_small+n_remote_small))
                allocate(rclo(n_local_small+n_remote_small))
                allocate(dclo(n_local_small+n_remote_small))

                isma = -1
                iclo = -1
                rclo = -1
                dclo = sum(extent ** 2)

                call near%set_lbound

                j = 0
                do n = 1, pcont%local_num + n_remote_small
                    call near%set_ubound(n)

                    if (pcont%is_small(n)) then
                        j = j + 1
                        isma(j) = n
                    endif
                enddo

                !---------------------------------------------------------------------
                ! Determine locally closest parcel:
                call find_closest_parcel_locally(pcont, n_local_small, n_remote_small, isma, iclo, rclo, dclo)
            endif

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(pcont, n_local_small, iclo, rclo, dclo)

            !---------------------------------------------------------------------
            ! We only need to check up to n_local_small because after the call to
            ! find_closest_parcel_globally all small parcels know their closest
            ! neighbour on the original MPI rank.
            if (.not. l_no_small) then
                do n = 1, n_local_small
                    if (iclo(n) == 0) then
                        write(*,*) iclo(n), dclo(n)
                        call mpi_exit_on_error('Merge error: no near enough neighbour found.')
                    endif
                enddo
            endif

            call stop_timer(merge_nearest_timer)

            !---------------------------------------------------------------------
            ! Figure out the mergers:
            if (subcomm%comm /= MPI_COMM_NULL) then
                call resolve_tree(isma, iclo, rclo, n_local_small)
            endif

            timings(merge_nearest_timer)%n_calls = timings(merge_nearest_timer)%n_calls - 1

            call start_timer(merge_nearest_timer)

            if (.not. l_no_small) then
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
            endif

            call deallocate_parcel_id_buffers

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: It cannot happen that the closest parcel is a small parcel
            !       on another MPI rank that is sent elsewhere.
            call gather_remote_parcels(pcont, n_local_small, n_invalid, rclo, iclo, isma, inva)

            !------------------------------------------------------------------
            ! Sanity check: Indices of close parcels must be smaller equal to the
            ! number of local parcels:
            if (.not. l_no_small) then
                if (maxval(iclo(1:n_local_small)) > pcont%local_num) then
                    call mpi_exit_on_error(&
                        "in parcel_nearest::find_nearest: Close parcel index out of range.")
                endif
            endif

            if (allocated(rclo)) then
                deallocate(rclo)
            endif

            if (allocated(dclo)) then
                deallocate(dclo)
            endif

            call near%dealloc

            call tree%finalise

            call stop_timer(merge_nearest_timer)

        end subroutine find_nearest

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! If a parcel communication happens across a periodic boundary, we must
        ! shift the parcel position such that it ends up in the sub-domain.
        subroutine apply_periodic_shift(pos, dir)
            double precision, intent(inout) :: pos(:)
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

        ! Find the nearest grid point to each small parcel (to be merged)
        ! and search over the surrounding 8 grid cells for the closest parcel.
        ! This operation is performed locally.
        subroutine find_closest_parcel_locally(pcont, n_local_small, n_remote_small, isma, iclo, rclo, dclo)
            class(pc_type),   intent(in)    :: pcont
            integer,          intent(in)    :: n_local_small
            integer,          intent(in)    :: n_remote_small
            integer,          intent(inout) :: isma(0:)
            integer,          intent(inout) :: iclo(:)
            integer,          intent(inout) :: rclo(:)
            double precision, intent(inout) :: dclo(:)
            integer                         :: ix0, iy0, iz0, ijk, n, ix, iy, iz, m, k, ic, is
            double precision                :: delx, dely, delz, dsq, dsqmin, x_small, y_small, z_small


            do m = 1, n_local_small + n_remote_small
                is = isma(m)

                x_small = pcont%position(1, is)
                y_small = pcont%position(2, is)
                z_small = pcont%position(3, is)
                ! Parcel "is" is small and should be merged; find closest other:
                ix0 = nint(dxi(1) * (x_small - box%halo_lower(1))) ! ranges from 0 to box%halo_size(1)
                iy0 = nint(dxi(2) * (y_small - box%halo_lower(2))) ! ranges from 0 to box%halo_size(2)
                iz0 = nint(dxi(3) * (z_small - box%lower(3)))      ! ranges from 0 to nz

                ! Grid point (ix0, iy0, iz0) is closest to parcel "is"

                dsqmin = sum(extent ** 2)
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
                            do k = near%kc1(ijk), near%kc2(ijk)
                                n = near%node(k)
                                if (n .ne. is) then
                                    delz = pcont%position(3, n) - z_small
                                    if (delz * delz < dsqmin) then
                                        delx = get_delx(pcont%position(1, n), x_small)
                                        dely = get_dely(pcont%position(2, n), y_small)
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

                ! Store the MPI rank and the index of the parcel to be potentially merged with:
                isma(m) = is
                iclo(m) = ic
                rclo(m) = cart%rank
                dclo(m) = dsqmin
                tree%l_merged(is) = .false.
                if (ic > 0) then
                    tree%l_merged(ic) = .false.
                endif
            enddo

            call update_remote_indices(pcont, n_local_small, n_remote_small, isma, iclo, rclo, dclo)


        end subroutine find_closest_parcel_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine update_remote_indices(pcont, n_local_small, n_remote_small, isma, iclo, rclo, dclo)
            class(pc_type),   intent(in)    :: pcont
            integer,          intent(in)    :: n_local_small
            integer,          intent(in)    :: n_remote_small
            integer,          intent(inout) :: isma(0:)
            integer,          intent(inout) :: iclo(:)
            integer,          intent(inout) :: rclo(:)
            double precision, intent(inout) :: dclo(:)
            integer                         :: m, is, ic

            !---------------------------------------------------------------------
            ! Update isma, iclo and rclo with indices of remote parcels:
            do m = 1, n_local_small + n_remote_small
                is = isma(m)
                ic = iclo(m)

                if ((is > pcont%local_num) .and. (ic > pcont%local_num)) then
                    ! A remote small parcel points to another remote small parcel. The remotes do not necessarily
                    ! need to be the same. Also, it can be a dual-link. As the same distance is evaluated on
                    ! the other two MPI ranks, *this* MPI rank must set the distance between the parcels to the
                    ! maximum value as otherwise the function "find_closest_parcel_globally" may think the
                    ! parcels belong to *this* MPI rank due to round-offs in the distance calculation. The
                    ! function "find_closest_parcel_globally" always sets "rclo" to the MPI source.
                    dclo(m) = huge(0.0d0) ! huge(x) returns the maximum value of this type
                else if (ic > pcont%local_num) then
                    ! A local small parcel points to a remote small parcel.
                    ! The index *ic* is larger than the local number of parcels
                    ! we must update *iclo(m)* and *rclo(m)* with the index of the parcel stored
                    ! on the other MPI rank.
                    iclo(m) = near%pidsma(ic)
                    rclo(m) = near%rsma(ic)
                endif
            enddo

            !------------------------------------------------------------------
            ! Sanity check: Indices of small parcels must be smaller equal to the
            ! number of local parcels:
            do m = 1, n_local_small
                if (isma(m) > pcont%local_num) then
                    call mpi_exit_on_error(&
                        'in in parcel_nearest::update_remote_indices: Small parcel index out of range.')
                endif
            enddo
        end subroutine update_remote_indices

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! For all small parcels in boundary cells we must now determine the
        ! closest parcel globally (that is we must compare with the neighbour rank).
        ! For this purpose we send the information of the small parcels that we received
        ! back to the original ranks which then figure out the closest parcel.
        ! Note: The information about the received small parcels is stored in the last
        !       n_remote_small entries of isma, iclo, rclo and dlco
        subroutine find_closest_parcel_globally(pcont, n_local_small, iclo, rclo, dclo)
            class(pc_type),   intent(in)            :: pcont
            integer,          intent(in)            :: n_local_small
            integer,          intent(inout)         :: iclo(:)
            integer,          intent(inout)         :: rclo(:)
            double precision, intent(inout)         :: dclo(:)
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: recv_size, send_size, buf_sizes(8)
            integer                                 :: tag, recv_count, n, l, i, m, k, j
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
                    k = j + pcont%local_num + l


                    ! merge index on *this* rank
                    m = n_local_small + j + l

                    send_buf(i)   = dble(near%midsma(k)) ! merge index on remote rank
                    send_buf(i+1) = dclo(m)              ! distance to closest parcel
                    send_buf(i+2) = dble(iclo(m))        ! parcel index of closest parcel
                enddo

                j = j + n_neighbour_small(n)

                ! receive order of small parcels
                tag = RECV_NEIGHBOUR_TAG(n)

#ifndef NDEBUG
                if (small_recv_order(n) /= tag) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::find_closest_parcel_globally: Wrong message order.")
                endif

                if (small_recv_count(n) /= n_neighbour_small(n)) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::find_closest_parcel_globally: Wrong number of parcels.")
                endif
#endif

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(tag)%rank,    &
                               SEND_NEIGHBOUR_TAG(tag), &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of parcel_nearest::find_closest_parcel_globally.")

            enddo

            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size,               &
                                           cart)

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::find_closest_parcel_globally: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              recv_status,              &
                              cart%err)

                call mpi_check_for_error(cart, &
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
                        rclo(m) = neighbours(n)%rank
                    endif
                enddo

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of parcel_nearest::find_closest_parcel_globally.")

            call deallocate_mpi_buffers

        end subroutine find_closest_parcel_globally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! https://github.com/mpi-forum/mpi-forum-historic/issues/413
        ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node294.htm
        ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node279.htm
        subroutine resolve_tree(isma, iclo, rclo, n_local_small)
            integer, intent(inout) :: isma(0:)
            integer, intent(inout) :: iclo(:)
            integer, intent(inout) :: rclo(:)
            integer, intent(inout) :: n_local_small
            integer                :: ic, rc, is, m, j
            logical                :: l_helper
            logical                :: l_continue_iteration, l_do_merge(n_local_small)
            logical                :: l_isolated_dual_link(n_local_small)

            call start_timer(merge_tree_resolve_timer)

            !------------------------------------------------------------------
            ! Exchange information:

            call tree%gather_info(iclo, rclo, n_local_small)

            !------------------------------------------------------------------
            ! Resolve tree now:

            ! First, iterative, stage
            l_continue_iteration = .true.

            do while (l_continue_iteration)
                l_continue_iteration = .false.
                ! reset relevant properties for candidate mergers

                do m = 1, n_local_small
                    is = isma(m)
                    ! only consider links that still may be merging
                    ! reset relevant properties
                    if (.not. tree%l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)
                        tree%l_leaf(is) = .true.
                        call tree%put_avail(rc, ic, .true.)
                    endif
                enddo

                ! Exchange information:
                call tree%sync_avail

                ! determine leaf parcels
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. tree%l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)
                        call tree%put_leaf(rc, ic, .false.)
                    endif
                enddo

                ! Exchange information:
                call tree%sync_leaf

                ! filter out parcels that are "unavailable" for merging
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. tree%l_merged(is)) then
                        if (.not. tree%l_leaf(is)) then
                            ic = iclo(m)
                            rc = rclo(m)
                            call tree%put_avail(rc, ic, .false.)
                        endif
                    endif
                enddo

                ! Exchange information:
                call tree%sync_avail

                ! identify mergers in this iteration
                do m = 1, n_local_small
                    is = isma(m)

                    if (.not. tree%l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        l_helper = tree%get_avail(rc, ic)

                        if (tree%l_leaf(is) .and. l_helper) then
                            l_continue_iteration = .true. ! merger means continue iteration
                            tree%l_merged(is) = .true.
                            call tree%put_merged(rc, ic, .true.)
                        endif
                    endif
                enddo

                ! Exchange information:
                call tree%sync_merged

                call start_timer(nearest_allreduce_timer)
                ! Performance improvement: We actually only need to synchronize with neighbours
                call MPI_Allreduce(MPI_IN_PLACE,            &
                                   l_continue_iteration,    &
                                   1,                       &
                                   MPI_LOGICAL,             &
                                   MPI_LOR,                 &
                                   subcomm%comm,            &
                                   subcomm%err)
                call stop_timer(nearest_allreduce_timer)
                call mpi_check_for_error(subcomm, &
                    "in MPI_Allreduce of parcel_nearest::resolve_tree.")
            enddo

            ! No barrier necessary because of the blocking MPI_Allreduce that acts like
            ! a barrier!

            ! Second stage, related to dual links
            do m = 1, n_local_small
                is = isma(m)

                if (.not. tree%l_merged(is)) then
                    if (tree%l_leaf(is)) then ! set in last iteration of stage 1
                        ic = iclo(m)
                        rc = rclo(m)
                        call tree%put_avail(rc, ic, .true.)
                    endif
                endif
            enddo

            ! Exchange information:
            call tree%sync_avail

            ! Second stage
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                l_do_merge(m) = .false.
                l_isolated_dual_link(m) = .false.

                if (tree%l_merged(is) .and. tree%l_leaf(is)) then
                    ! previously identified mergers: keep
                    l_do_merge(m) = .true.
                    !----------------------------------------------------------
                    ! begin of sanity check
                    ! After first stage mergers parcel cannot be both initiator
                    ! and receiver in stage 1
                    l_helper = tree%get_leaf(rc, ic)

                    if (l_helper) then
                        call mpi_exit_on_error(&
                            'in parcel_nearest::resolve_tree: First stage error')
                    endif

                    ! end of sanity check
                    !----------------------------------------------------------

                elseif (.not. tree%l_merged(is)) then
                    if (tree%l_leaf(is)) then
                        ! links from leafs
                        l_do_merge(m) = .true.
                    elseif (.not. tree%l_available(is)) then
                        ! Above means parcels that have been made 'available' do not keep outgoing links

                        l_helper = tree%get_avail(rc, ic)

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
                            if (cart%rank <= rc) then
                                ! The MPI rank with lower number makes its parcel
                                ! available.
                                tree%l_available(is) = .true.
                            endif
                        endif
                    endif
                endif
            enddo

            ! Exchange information:
            call tree%sync_avail

            !------------------------------------------------------
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)

                if ((l_do_merge(m) .eqv. .false.) .and. l_isolated_dual_link(m)) then
                    ! isolated dual link

                    l_helper = tree%get_avail(rc, ic)

                    if (l_helper) then
                        ! merge this parcel into ic along with the leaf parcels
                        l_do_merge(m) = .true.
                    !else
                    !   ! Dual link is resolved on other rank
                    endif
                endif
                !------------------------------------------------------
            enddo

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

            call tree%free_memory

            call stop_timer(merge_tree_resolve_timer)
        end subroutine resolve_tree

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This routine fills the *pid* buffers with 2 entries per small parcel.
        ! The first entry is the local parcel index in the parcel container and
        ! the second entry is the merge index (usually accessed with *m*).
        subroutine locate_parcel_in_boundary_cell(pcont, m, n)
            class(pc_type), intent(in) :: pcont
            integer,        intent(in) :: m, n
            integer                    :: k, ix, iy

            ! nearest global grid point
            ix = nint(dxi(1) * (pcont%position(1, n) - lower(1)))
            iy = nint(dxi(2) * (pcont%position(2, n) - lower(2)))

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
        subroutine send_small_parcel_bndry_info(pcont, n_remote_small)
            class(pc_type),           intent(inout) :: pcont
            integer,                  intent(inout) :: n_remote_small
            integer,          dimension(:), pointer :: send_ptr
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: recv_size, send_size
            integer                                 :: recv_count, n, i ,j, l, m, pid, k
            integer, parameter                      :: n_entries = 5
            integer, allocatable                    :: rtmp(:), pidtmp(:), midtmp(:)
            ! rtmp: MPI rank remote small parcel belongs to
            ! pidtmp: parcel index of remote parcel (on the owning rank)
            ! midtmp: m index of remote parcel (on the owning rank)

            n_neighbour_small = 0

            ! dummy allocate; will be updated on the fly
            allocate(near%pidsma(0))
            allocate(near%rsma(0))
            allocate(near%midsma(0))

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
                        send_buf(i:i+2) = pcont%position(:, pid)
                        send_buf(i+3) = dble(pid)
                        send_buf(i+4) = dble(m)
                    enddo
                endif

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of parcel_nearest::send_small_parcel_bndry_info.")
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size,               &
                                           cart)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              recv_status,              &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Recv of parcel_nearest::send_small_parcel_bndry_info.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::send_small_parcel_bndry_info: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

#ifndef NDEBUG
                small_recv_order(n) = RECV_NEIGHBOUR_TAG(n)
                small_recv_count(n) = recv_count
#endif

                i = pcont%local_num+1
                j = pcont%local_num+1+size(near%rsma) + recv_count
                allocate(rtmp(i:j))
                allocate(pidtmp(i:j))
                allocate(midtmp(i:j))

                ! copy old over
                if (size(near%rsma) > 0) then
                    rtmp(pcont%local_num+1:pcont%local_num+size(near%rsma)) = near%rsma
                    pidtmp(pcont%local_num+1:pcont%local_num+size(near%pidsma)) = near%pidsma
                    midtmp(pcont%local_num+1:pcont%local_num+size(near%midsma)) = near%midsma
                endif

                deallocate(near%pidsma)
                deallocate(near%rsma)
                deallocate(near%midsma)

                if (recv_count > 0) then
                    ! unpack parcel position and parcel index to recv buffer
                    do l = 1, recv_count
                        i = 1 + (l-1) * n_entries
                        k = sum(n_neighbour_small) + pcont%local_num + l
                        pcont%position(:, k) = recv_buf(i:i+2)
                        pcont%volume(k) = zero    ! set volume / area to zero as each parcel is small
                        pidtmp(k) = nint(recv_buf(i+3))
                        rtmp(k) = neighbours(n)%rank
                        midtmp(k) = nint(recv_buf(i+4))

                        !------------------------------------------------------
                        ! We must also assign incoming small parcels to cells
                        ! Note: We must apply a shift for parcels communicated
                        !       across a periodic boundary.
                        call apply_periodic_shift(pcont%position(:, k), &
                                                  RECV_NEIGHBOUR_TAG(n))

                        ! Now we can safely assign the local cell index:
                        call near%parcel_to_local_cell_index(pcont, k)
                    enddo
                    n_neighbour_small(n) = n_neighbour_small(n) + recv_count
                endif

                ! *move_alloc* deallocates pidtmp, rtmp and midtmp
                call move_alloc(pidtmp, near%pidsma)
                call move_alloc(rtmp, near%rsma)
                call move_alloc(midtmp, near%midsma)

                deallocate(recv_buf)
            enddo

            n_remote_small = sum(n_neighbour_small)

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of parcel_nearest::send_small_parcel_bndry_info.")


            call deallocate_parcel_buffers

        end subroutine send_small_parcel_bndry_info

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine gather_remote_parcels(pcont, n_local_small, n_invalid, rclo, iclo, isma, inva)
            class(pc_type),              intent(inout) :: pcont
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
            integer                                    :: recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(pcont%attr_num)
            integer                                    :: m, rc, ic, is, n, i, j, k, iv
            integer                                    :: n_entries
            integer                                    :: n_registered(8)
#ifndef NDEBUG
            integer(kind=int64)                        :: n_total
#endif

            !------------------------------------------------------------------
            ! Figure out how many parcels we send and allocate buffers:
            n_parcel_sends = 0

            do m = 1, n_local_small
                rc = rclo(m)
                if (cart%rank /= rc) then
                    n = get_neighbour_from_rank(rc)
                    n_parcel_sends(n) = n_parcel_sends(n) + 1
                endif
            enddo

            n_invalid = sum(n_parcel_sends)

            ! We must send all parcel attributes (attr_num) plus
            ! the index of the close parcel ic (1)
            n_entries = pcont%attr_num + 1
            n_registered = n_parcel_sends * n_entries
            call allocate_mpi_buffers(n_registered)

            !------------------------------------------------------------------
            ! Fill all buffers:
            n_registered = 0
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                inva(m) = is ! copy initial isma to inva for parcel deletion
                if (cart%rank /= rc) then
                    ! The closest parcel to this small parcel *is*
                    ! is on another MPI rank. We must send this parcel
                    ! to that rank.

                    n = get_neighbour_from_rank(rc)

                    call get_mpi_buffer(n, send_buf)

                    n_registered(n) = n_registered(n) + 1

                    call pcont%serialize(is, buffer)
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

            iv = 0
            if (n_local_small > 0) then
                ! We must now remove all invalid entries in isma and
                ! iclo and also update the value of n_local_small:
                iv = n_local_small  ! we need to keep the original value for filling *inva*
                n_local_small = count(isma(1:iv) /= -1)
                isma(1:n_local_small) = pack(isma(1:iv), isma(1:iv) /= -1)
                iclo(1:n_local_small) = pack(iclo(1:iv), iclo(1:iv) /= -1)
                rclo(1:n_local_small) = pack(rclo(1:iv), rclo(1:iv) /= -1)
            endif

            !------------------------------------------------------------------
            ! Communicate parcels:
            do n = 1, 8
                call get_mpi_buffer(n, send_buf)

                send_size = n_parcel_sends(n) * n_entries

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of parcel_nearest::gather_remote_parcels.")
            enddo

            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size,               &
                                           cart)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              recv_status,              &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Recv of parcel_nearest::gather_remote_parcels.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::gather_remote_parcels: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                if (recv_count > 0) then
                    ! Set the current value of the *m* index to
                    ! its current maximum value (excluding small
                    ! parcel received earlier from other MPI ranks.
                    m = n_local_small
                    do k = 1, recv_count
                        ! We receive a small parcel;
                        ! append it to the container with index
                        ! "pcont%local_num+1"
                        pcont%local_num = pcont%local_num + 1
                        j = n_entries * k
                        i = j - n_entries + 1
                        buffer = recv_buf(i:j-1)
                        call pcont%deserialize(pcont%local_num, buffer)

                        ! Add the small parcel to isma and inva, and
                        ! its closest parcel to iclo
                        m = m + 1
                        iclo(m) = nint(recv_buf(j))
                        isma(m) = pcont%local_num
                        iv = iv + 1
                        inva(iv) = pcont%local_num
                    enddo
                    ! The last value of *m* is our new number of
                    ! small parcels. However, this still includes
                    ! invalid entries in isma and iclo that we will
                    ! remove next.
                    n_local_small = m
                endif

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of parcel_nearest::gather_remote_parcels.")

            call deallocate_mpi_buffers

#ifndef NDEBUG
            if (subcomm%comm == world%comm) then
                n_total = pcont%local_num - n_invalid
                call mpi_blocking_reduce(n_total, MPI_SUM, world)
                if ((world%rank == world%root) .and. (.not. n_total == pcont%total_num)) then
                    call mpi_exit_on_error(&
                        "in parcel_nearest::gather_remote_parcels: We lost parcels.")
                endif
            endif
#endif
        end subroutine gather_remote_parcels

end module parcel_nearest
