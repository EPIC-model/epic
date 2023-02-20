!==============================================================================
!               Finds the parcels nearest every "small" parcel
!==============================================================================
module parcel_nearest
    use constants, only : zero !pi, f12
    use parcel_container, only : parcels, n_parcels, get_delx, get_dely, n_par_attrib
    use parameters, only : dx, dxi, vcell, hli, lower, extent, ncell, nx, ny, nz, vmin, max_num_parcels
    use options, only : parcel
    use timer, only : start_timer, stop_timer
    use mpi_communicator
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_error, mpi_check_for_message
    use iso_c_binding, only : c_sizeof
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
                         , get_parcel_buffer_ptr        &
                         , communicate_parcels
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
    logical, allocatable :: l_merged(:)    ! indicates parcels merged in first stage

#ifndef NDEBUG
    ! Logicals that are only needed for sanity checks
    logical, allocatable :: l_is_merged(:) ! SANITY CHECK ONLY
    logical, allocatable :: l_small(:)     ! SANITY CHECK ONLY
    logical, allocatable :: l_close(:)     ! SANITY CHECK ONLY
#endif

    integer              :: n_neighbour_small(8)  ! number of small parcels received
    integer              :: n_remote_small        ! sum(n_neighbour_small)
    integer, allocatable :: remote_small_rank(:)  ! rank of small parcel received
    integer, allocatable :: remote_small_pid(:)   ! index of remote small parcel
    integer, allocatable :: remote_m_index(:)     ! *m* in *isma(m)*

    logical :: l_continue_iteration, l_do_merge

    type(MPI_Win) :: win_merged, win_avail, win_leaf

#ifndef NDEBUG
    type(MPI_Win) :: win_is_merged, win_small, win_close
#endif

    public :: find_nearest, merge_nearest_timer, merge_tree_resolve_timer

    contains

        subroutine nearest_allocate
            integer (KIND=MPI_ADDRESS_KIND) :: win_size
            integer, parameter              :: disp_unit = c_sizeof(l_do_merge) ! size of logical in bytes

            if (.not. allocated(nppc)) then
                allocate(nppc(box%halo_ncell))
                allocate(kc1(box%halo_ncell))
                allocate(kc2(box%halo_ncell))
                allocate(loca(max_num_parcels))
                allocate(node(max_num_parcels))
                allocate(l_leaf(max_num_parcels))
                allocate(l_available(max_num_parcels))
                allocate(l_merged(max_num_parcels))

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
                call MPI_Win_create(l_leaf,         &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_leaf,       &
                                    comm%err)
                call MPI_Win_create(l_available,    &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_avail,      &
                                    comm%err)
                call MPI_Win_create(l_merged,       &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_merged,     &
                                    comm%err)

#ifndef NDEBUG
                allocate(l_is_merged(max_num_parcels))
                allocate(l_small(max_num_parcels))
                allocate(l_close(max_num_parcels))

                call MPI_Win_create(l_is_merged,    &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_is_merged,  &
                                    comm%err)
                call MPI_Win_create(l_small,        &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_small,      &
                                    comm%err)
                call MPI_Win_create(l_close,        &
                                    win_size,       &
                                    disp_unit,      &
                                    MPI_INFO_NULL,  &
                                    comm%world,     &
                                    win_close,      &
                                    comm%err)
#endif
            endif
        end subroutine nearest_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine nearest_deallocate
            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
                deallocate(l_leaf)
                deallocate(l_available)
                deallocate(l_merged)
#ifndef NDEBUG
                deallocate(l_is_merged)
                deallocate(l_small)
                deallocate(l_close)
#endif
            endif
        end subroutine nearest_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! @param[out] isma indices of small parcels
        ! @param[out] iclo indices of close parcels
        ! @param[out] n_local_small the array size of isma and iclo
        ! @post
        !   - isma must be sorted in ascending order
        !   - isma and iclo must be filled contiguously
        !   - parcel indices in isma cannot be in iclo, and vice-versa
        !   - the m-th entry in isma relates to the m-th entry in iclo
        subroutine find_nearest(isma, iclo, n_local_small)
            integer, allocatable, intent(out) :: isma(:)
            integer, allocatable, intent(out) :: iclo(:)
            integer,              intent(out) :: n_local_small
            integer                           :: n_global_small
            integer                           :: ijk, n, ix, iy, k, j
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
            n_local_small = 0
            n_global_small = 0
            n_neighbour_small = 0
            n_remote_small = 0
            nppc = 0 !nppc(ijk) will contain the number of parcels in grid cell ijk

            ! Bin parcels in cells:
            ! Form list of small parcels:
            do n = 1, n_parcels

                call parcel_to_local_cell_index(n, ix, iy)

                if (parcels%volume(n) < vmin) then
                    n_local_small = n_local_small + 1

                    ! If a small parcel is in a boundary cell, a duplicate must
                    ! be sent to the neighbour rank. This call checks if the parcel
                    ! must be sent and fills the send buffers.
                    call locate_parcel_in_boundary_cell(n_local_small, n, ix, iy)
                endif
            enddo

            call MPI_Allreduce(n_local_small, n_global_small, 1, MPI_INTEGER, MPI_SUM, comm%world, comm%err)

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
            call send_small_parcel_bndry_info(n_local_small)

            ! We must also assign incoming small parcels to cells
            ! Note: We must apply a shift for parcels communicated
            !       across a periodic boundary.
            do n = n_parcels + 1, n_parcels + n_remote_small
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
            allocate(iclo(n_local_small + n_remote_small))
            allocate(rclo(n_local_small + n_remote_small))
            allocate(dclo(n_local_small + n_remote_small))

            isma = 0
            iclo = 0
            rclo = comm%rank
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
#ifndef NDEBUG
                l_is_merged(n) = .false.! SANITY CHECK ONLY
                l_small(n) = .false. ! SANITY CHECK ONLY
                l_close(n) = .false. ! SANITY CHECK ONLY
#endif
            enddo

            !---------------------------------------------------------------------
            ! Determine locally closest parcel:
            call find_closest_parcel_locally(n_local_small, isma, iclo, rclo, dclo)

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(n_local_small, iclo, rclo, dclo)

#ifndef NDEBUG
            write(*,*) 'start merging, n_local_small='
            write(*,*) n_local_small
#endif

            call stop_timer(merge_nearest_timer)

            !---------------------------------------------------------------------
            ! Figure out the mergers:
            call resolve_tree(isma, iclo, rclo, n_local_small)

            call deallocate_parcel_id_buffers

            ! We must store the parcel index of each small parcel.
            ! We need to allocate the invalid buffer, therefore
            ! the second argument is .true.
            call allocate_parcel_id_buffers(1, .true.)

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: If the closest parcel is a small parcel on another MPI rank,
            !       we must check if this remote parcel is not sent elsewhere.
            !       FIXME Double-check: This might actually never happen.
            call gather_remote_parcels(n_local_small, rclo, iclo, isma)

            call deallocate_parcel_id_buffers

            call nearest_deallocate

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
                dclo(m) = dsqmin
                l_merged(is) = .false.
                l_merged(ic) = .false.
            enddo

            ! Update isma, iclo and rclo with indices of remote parcels:
            do m = n_local_small + 1, n_local_small + n_remote_small
                is = isma(m)
                ic = iclo(m)
                ! If the index *is* or *ic* are larger than the local number of parcels
                ! we must update *isma(m)* or *iclo(m)* with the index of the parcel stored
                ! on the other MPI rank. Otherwise keep the current value. This operation is
                ! achieved with the built-in Fortran routine *merge*.
                isma(m) = merge(remote_small_pid(m), is, is > n_parcels)
                iclo(m) = merge(remote_small_pid(m), ic, ic > n_parcels)
                rclo(m) = merge(remote_small_rank(m), comm%rank, ic > n_parcels)
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
            type(MPI_Request)                       :: request
            type(MPI_Status)                        :: recv_status, send_status
            integer                                 :: recv_size, send_size
            integer                                 :: tag, source, recv_count, n, l, i, m
            integer, parameter                      :: n_entries = 3


            !------------------------------------------------------------------
            ! Communicate with neighbours:
            do n = 1, 8

                ! we send the distance, the remote parcel index and the
                ! remote merge index *m* --> n_entries = 3.
                send_size = n_neighbour_small(n) * n_entries

                allocate(send_buf(send_size))

                do l = 1, n_neighbour_small(n)
                    i = 1 + (l-1) * n_entries
                    m = n_local_small + 1 + sum(n_neighbour_small(1:n-1))
                    send_buf(i)   = remote_m_index(m)       ! merge_index
                    send_buf(i+1) = dclo(m)                 ! distance to closest parcel
                    send_buf(i+2) = dble(iclo(m))           ! parcel index of closest parcel
                enddo

                call MPI_Isend(send_buf,                &
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

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_nearest::find_closest_parcel_globally: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf,                 &
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

                deallocate(send_buf)
                deallocate(recv_buf)

                call mpi_check_for_error(&
                    "in MPI_Wait of parcel_nearest::find_closest_parcel_globally.")
            enddo

        end subroutine find_closest_parcel_globally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine resolve_tree(isma, iclo, rclo, n_local_small)
            integer, intent(inout)         :: isma(0:)
            integer, intent(inout)         :: iclo(:)
            integer, intent(inout)         :: rclo(:)
            integer, intent(inout)         :: n_local_small
            integer                        :: ic, rc, is, m, j
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

            do while (l_continue_iteration)
                l_continue_iteration = .false.
                ! reset relevant properties for candidate mergers

                call MPI_Win_fence(0, win_avail, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

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
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         icm,               &
                                         max_num_parcels,   &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)
                            call mpi_check_for_error(&
                                "in MPI_Put of parcel_nearest::resolve_tree.")
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_avail, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

                call MPI_Win_fence(0, win_leaf, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

                ! determine leaf parcels
                do m = 1, n_local_small
                    is = isma(m)
                    if (.not. l_merged(is)) then
                        ic = iclo(m)
                        rc = rclo(m)

                        if (rc == comm%rank) then
                            l_leaf(ic) = .false.
                        else
                            l_helper = .false.
                            icm = ic
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         icm,               &
                                         max_num_parcels,   &
                                         MPI_LOGICAL,       &
                                         win_leaf,          &
                                         comm%err)
                            call mpi_check_for_error(&
                                "in MPI_Put of parcel_nearest::resolve_tree.")
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_leaf, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

                call MPI_Win_fence(0, win_avail, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

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
                                l_helper = .false.
                                icm = ic
                                call MPI_Put(l_helper,          &
                                             1,                 &
                                             MPI_LOGICAL,       &
                                             rc,                &
                                             icm,               &
                                             max_num_parcels,   &
                                             MPI_LOGICAL,       &
                                             win_avail,         &
                                             comm%err)
                                call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")
                            endif
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_avail, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

                call MPI_Win_fence(0, win_merged, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

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
                            !     MPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                            !         target_disp, target_count, target_datatype, win, ierror)
                            !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: origin_addr
                            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
                            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
                            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
                            !     TYPE(MPI_Win), INTENT(IN) :: win
                            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
                            icm = ic
                            call MPI_Get(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         icm,               &
                                         max_num_parcels,   &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)
                            call mpi_check_for_error(&
                                    "in MPI_Get of parcel_nearest::resolve_tree.")
                        endif

                        if (l_leaf(is) .and. l_helper) then
                            l_continue_iteration = .true. ! merger means continue iteration
                            l_merged(is) = .true.

                            if (rc == comm%rank) then
                                l_merged(ic) = .true.
                            else
                                l_helper = .true.
                                icm = ic
                                call MPI_Put(l_helper,          &
                                            1,                  &
                                            MPI_LOGICAL,        &
                                            rc,                 &
                                            icm,                &
                                            max_num_parcels,    &
                                            MPI_LOGICAL,        &
                                            win_merged,         &
                                            comm%err)
                                call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")
                            endif
                        endif
                    endif
                enddo

                call MPI_Win_fence(0, win_merged, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

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

            call MPI_Win_fence(0, win_avail, comm%err)

            call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

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
                            l_helper = .true.
                            icm = ic
                            call MPI_Put(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         icm,               &
                                         max_num_parcels,   &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)

                            call mpi_check_for_error(&
                                    "in MPI_Put of parcel_nearest::resolve_tree.")
                        endif
                    endif
                endif
            enddo

            call MPI_Win_fence(0, win_avail, comm%err)

            call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

            call MPI_Win_fence(0, win_avail, comm%err)

            call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

            ! Second stage (hard to parallelise with openmp?)
            j = 0
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                l_do_merge = .false.
                if (l_merged(is) .and. l_leaf(is)) then
                    ! previously identified mergers: keep
                    l_do_merge=.true.
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

                        if (rc == comm%rank) then
                            l_helper = l_available(ic)
                        else
                            call MPI_Get(l_helper,          &
                                         1,                 &
                                         MPI_LOGICAL,       &
                                         rc,                &
                                         icm,               &
                                         max_num_parcels,   &
                                         MPI_LOGICAL,       &
                                         win_avail,         &
                                         comm%err)

                            call mpi_check_for_error(&
                                "in MPI_Get of parcel_nearest::resolve_tree.")
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
                            if (comm%rank <= rc) then
                                ! The MPI rank with lower number makes its parcel
                                ! available.
                                l_available(is) = .true.
                            endif
                        endif
                    endif
                endif

                call MPI_Win_fence(0, win_avail, comm%err)

                call mpi_check_for_error(&
                    "in MPI_Win_fence of parcel_nearest::resolve_tree.")

                if (l_do_merge) then
                    j = j + 1
                    isma(j) = is
                    iclo(j) = ic
                    rclo(j) = rc
#ifndef NDEBUG
                    l_is_merged(is) = .true.
                    l_is_merged(ic) = .true.
                    l_small(is) = .true.
                    l_close(ic) = .true.
#endif
                endif
            enddo
            n_local_small = j

#ifndef NDEBUG
            write(*,*) 'after second stage, n_local_small='
            write(*,*) n_local_small
            write(*,*) 'finished'

            ! MORE SANITY CHECKS
            ! CHECK ISMA ORDER
            do m = 1, n_local_small
                if (.not. (isma(m) > isma(m-1))) then
                    write(*,*) 'isma order broken'
                endif
            enddo

            ! 1. CHECK RESULTING MERGERS
            do m = 1, n_local_small
                if (.not. l_is_merged(isma(m))) write(*,*) 'merge_error: isma(m) not merged, m=', m
                if (.not. l_is_merged(iclo(m))) write(*,*) 'merge_error: iclo(m) not merged, m=', m
                if (.not. l_small(isma(m))) write(*,*) 'merge_error: isma(m) not marked as small, m=', m
                if (.not. l_close(iclo(m))) write(*,*) 'merge_error: iclo(m) not marked as close, m=', m
                if (l_close(isma(m))) write(*,*) 'merge_error: isma(m) both small and close, m=', m
                if (l_small(iclo(m))) write(*,*) 'merge_error: iclo(m) both small and close, m=', m
            enddo

            ! 2. CHECK MERGING PARCELS
            do n = 1, n_parcels
                if (parcels%volume(n) < vmin) then
                    if (.not. l_is_merged(n)) then
                        write(*,*) 'merge_error: parcel n not merged (should be), n=', n
                    endif
                    if (.not. (l_small(n) .or. l_close(n))) then
                        write(*,*) 'merge_error: parcel n not small or close (should be), n=', n
                    endif
                    if (l_small(n) .and. l_close(n)) then
                        write(*,*) 'merge_error: parcel n both small and close, n=', n
                    endif
                else
                    if (l_small(n)) then
                        write(*,*) 'merge_error: parcel n small (should not be), n=', n
                    endif
                    if (l_is_merged(n) .and. (.not. l_close(n))) then
                        write(*,*) 'merge_error: parcel n merged (should not be), n=', n
                    endif
                endif
            enddo
#endif
            call stop_timer(merge_tree_resolve_timer)
        end subroutine resolve_tree

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This routine fills the *pid* buffers with 2 entries per small parcel.
        ! The first entry is the local parcel index in the parcel container and
        ! the second entry is the merge index (usually accessed with *m*).
        subroutine locate_parcel_in_boundary_cell(m, n, ix, iy)
            integer, intent(in) :: m, n, ix, iy
            integer             :: k

            ! check lower x-direction
            if (ix == box%lo(1)) then

                if (iy == box%lo(2)) then
                    ! parcel in southwest corner with
                    ! neighbours: west, south and southwest
                    k = n_parcel_sends(MPI_SOUTHWEST) + 1
                    southwest_pid(2*k-1) = n
                    southwest_pid(2*k) = m
                    n_parcel_sends(MPI_SOUTHWEST) = k

                else if (iy == box%hi(2)) then
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
            if (ix >= box%hi(1)) then

                if (iy == box%lo(2)) then
                    ! parcel in southeast corner with
                    ! neighbours: east, south and southeast
                    k = n_parcel_sends(MPI_SOUTHEAST) + 1
                    southeast_pid(2*k-1) = n
                    southeast_pid(2*k) = m
                    n_parcel_sends(MPI_SOUTHEAST) = k

                else if (iy == box%hi(2)) then
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

            if (iy >= box%hi(2)) then
                k = n_parcel_sends(MPI_NORTH) + 1
                north_pid(2*k-1) = n
                north_pid(2*k) = m
                n_parcel_sends(MPI_NORTH) = k
            endif
        end subroutine locate_parcel_in_boundary_cell

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Send position of small parcels in boundary region, their local index
        ! in the parcel container and their merge index (*m*) to neighbours
        subroutine send_small_parcel_bndry_info(n_local_small)
            integer ,                    intent(in) :: n_local_small
            integer,          dimension(:), pointer :: send_ptr
            double precision, dimension(:), pointer :: send_buf
            double precision, allocatable           :: recv_buf(:)
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: recv_size, send_size
            integer                                 :: tag, source, recv_count, n, i ,j, l, m, pid
            integer, parameter                      :: n_entries = 5
            integer, allocatable                    :: tmp_rank(:), tmp_pid(:), tmp_m_index(:)

            n_neighbour_small = 0

            ! dummy allocate; will be updated on the fly
            allocate(remote_small_pid(0))
            allocate(remote_small_rank(0))
            allocate(remote_m_index(0))

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
                         send_buf(i:i+3) = parcels%position(:, pid)
                         send_buf(i+4) = dble(pid)
                         send_buf(i+5) = dble(m)
                     enddo
                endif

                call MPI_Isend(send_buf,                &
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

                call MPI_Recv(recv_buf,                 &
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

                i = n_local_small + 1
                j = n_local_small + size(remote_small_rank) + recv_count
                allocate(tmp_rank(i:j))
                allocate(tmp_pid(i:j))
                allocate(tmp_m_index(i:j))
                deallocate(remote_small_pid)
                deallocate(remote_small_rank)
                deallocate(remote_m_index)

                if (recv_count > 0) then
                    ! unpack parcel position and parcel index to recv buffer
                    do l = 1, recv_count
                        i = 1 + (l-1) * n_entries
                        j = sum(n_neighbour_small) + n_parcels + l
                        parcels%position(:, j) = recv_buf(i:i+3)
                        tmp_pid(j) = nint(recv_buf(i+4))
                        tmp_rank(j) = source
                        tmp_m_index(j) = nint(recv_buf(i+5))
                    enddo
                    n_neighbour_small(source) = n_neighbour_small(source) + recv_count
                endif

                ! *move_alloc* deallocates tmp_pid and tmp_rank
                call move_alloc(tmp_pid, remote_small_pid)
                call move_alloc(tmp_rank, remote_small_rank)
                call move_alloc(tmp_m_index, remote_m_index)

                deallocate(recv_buf)
            enddo

            n_remote_small = sum(n_neighbour_small)

            allocate(remote_small_pid(n_local_small+1:n_local_small+n_remote_small))
            allocate(remote_small_rank(n_local_small+1:n_local_small+n_remote_small))
            allocate(remote_m_index(n_local_small+1:n_local_small+n_remote_small))

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

            call mpi_check_for_error(&
                "in MPI_Waitall of parcel_nearest::send_small_parcel_bndry_info.")


            call deallocate_parcel_buffers

        end subroutine send_small_parcel_bndry_info

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine gather_remote_parcels(n_local_small, rclo, iclo, isma)
            integer,                     intent(inout) :: n_local_small
            integer,                     intent(in)    :: rclo(:)       ! MPI rank of closest parcel
            integer,                     intent(inout) :: iclo(:)
            integer,                     intent(inout) :: isma(0:)
            integer,          dimension(:), pointer    :: send_pid
            double precision, dimension(:), pointer    :: send_buf
            double precision, allocatable              :: recv_buf(:)
            type(MPI_Request)                          :: requests(8)
            type(MPI_Status)                           :: recv_status, send_statuses(8)
            integer                                    :: tag, source, recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(n_par_attrib)
            integer                                    :: m, rc, ic, is, n, iv, i, j, k
            integer, parameter                         :: n_entries = n_par_attrib + 1

            ! We must send all parcel attributes (n_par_attrib) plus
            ! the index of the close parcel ic (1)
            call allocate_parcel_buffers(n_entries)

            n_parcel_sends = 0

            !------------------------------------------------------------------
            ! Figure out which small parcels we need to send:
            iv = 1
            do m = 1, n_local_small
                rc = rclo(m)
                ic = iclo(m)
                is = isma(m)

                if (.not. rc == comm%rank) then
                    ! The closest parcel to this smalll parcel *is*
                    ! is on another MPI rank. We must send this parcel
                    ! to that rank.
                    call get_index_periodic(parcels%position(:, ic), i, j, k)
                    n = get_neighbour(i, j)

                    call get_parcel_buffer_ptr(n, send_pid, send_buf)
                    n_parcel_sends(n) = n_parcel_sends(n) + 1

                    call parcel_serialize(is, buffer)
                    j = n_entries * iv
                    i = j - n_entries + 1
                    send_buf(i:j-1) = buffer
                    send_buf(j) = dble(ic)

                    n = n_parcel_sends(n)
                    send_pid(n) = is
                    invalid(iv) = is
                    iv = iv + 1

                    ! Mark this as invalid
                    isma(m) = -1
                    iclo(m) = -1
                endif
            enddo

            !------------------------------------------------------------------
            ! Communicate parcels:
            do n = 1, 8
                call get_parcel_buffer_ptr(n, send_pid, send_buf)

                send_size = n_parcel_sends(n) * n_entries

                call MPI_Isend(send_buf,                &
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

                call MPI_Recv(recv_buf,                 &
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

                    ! Add the small parcel to isma and
                    ! its closest parcel to iclo
                    ic = nint(recv_buf(j))
                    m = m + 1
                    iclo(m) = ic
                    isma(m) = n_parcels
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

            ! delete parcel that we sent
            n = sum(n_parcel_sends)
            call parcel_delete(invalid, n)

            call deallocate_parcel_buffers

            ! We must now remove all invalid entries in isma and
            ! iclo and also update the value of n_local_small:
            n_local_small = count(isma(1:) /= -1)
            isma(1:n_local_small) = pack(isma(1:), isma(1:) /= -1)
            iclo(1:n_local_small) = pack(iclo, iclo /= -1)

#ifndef NDEBUG
            n = n_parcels
            call mpi_blocking_reduce(n, MPI_SUM)
            if ((comm%rank == comm%master) .and. (.not. n == n_total_parcels)) then
                call mpi_exit_on_error(&
                    "in parcel_nearest::gather_remote_parcels: We lost parcels.")
            endif
#endif
        end subroutine gather_remote_parcels

end module parcel_nearest
