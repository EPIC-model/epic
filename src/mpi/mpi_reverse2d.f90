module mpi_reverse2d
    use mpi_environment
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_error
    implicit none

    private

    type :: reorder_type
        integer                       :: dir                            ! reorder index
        integer, allocatable          :: dest(:)                        ! destination rank
        integer, allocatable          :: send_recv_count(:)
        integer, allocatable          :: recv_count(:)                  ! number of receives of
                                                                        ! the reversing dimension
        integer, allocatable          :: send_offset(:), recv_offset(:)
        double precision, allocatable :: send_buffer(:), recv_buffer(:)

        ! arrays for halo fill
        double precision, allocatable :: lo_buffer(:), hi_halo_buffer(:)
        double precision, allocatable :: hi_buffer(:), lo_halo_buffer(:)
        integer                       :: lo_rank, hi_rank

    end type reorder_type

    logical :: l_initialised = .false.

    type(communicator) :: x_comm
    type(communicator) :: y_comm

    type(reorder_type) :: x_reo
    type(reorder_type) :: y_reo

    public :: reverse2d_x              &
            , reverse2d_y              &
            , initialise_mpi_reverse2d &
            , finalise_mpi_reverse2d

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine initialise_mpi_reverse2d
            if (l_initialised) then
                return
            endif
            l_initialised = .true.

            call initialise_reversing(x_reo, x_comm, 1)

            call initialise_reversing(y_reo, y_comm, 2)

        end subroutine initialise_mpi_reverse2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_mpi_reverse2d
            if (.not. l_initialised) then
                return
            endif

            call finalise_reversing(x_reo)

            call finalise_reversing(y_reo)

        end subroutine finalise_mpi_reverse2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_reversing(reo)
            type(reorder_type), intent(inout) :: reo

            deallocate(reo%send_recv_count)
            deallocate(reo%send_offset)
            deallocate(reo%recv_count)
            deallocate(reo%recv_offset)

            deallocate(reo%dest)

            deallocate(reo%send_buffer)
            deallocate(reo%recv_buffer)


            deallocate(reo%lo_buffer)
            deallocate(reo%hi_buffer)
            deallocate(reo%hi_halo_buffer)
            deallocate(reo%lo_halo_buffer)

        end subroutine finalise_reversing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine initialise_reversing(reo, sub_comm, dir)
            type(reorder_type), intent(inout) :: reo
            type(communicator), intent(inout) :: sub_comm
            integer,            intent(in)    :: dir
            integer, allocatable              :: rlos(:), rhis(:)
            integer                           :: i, j, lo, hi, rank, d
            logical                           :: remain_dims(2)

            reo%dir = dir

            ! if dir = 1 --> d = 2, i.e. if dir = x --> d = y
            ! if dir = 2 --> d = 1, i.e. if dir = y --> d = y
            d = mod(dir, 2) + 1

            remain_dims(dir) = .true.
            remain_dims(d) = .false.

            call MPI_Cart_sub(cart%comm, remain_dims, sub_comm%comm, cart%err)

            call MPI_Comm_rank(sub_comm%comm, sub_comm%rank, sub_comm%err)

            call MPI_Comm_size(sub_comm%comm, sub_comm%size, sub_comm%err)

            allocate(reo%send_recv_count(sub_comm%size))
            allocate(reo%send_offset(sub_comm%size))
            allocate(reo%recv_count(sub_comm%size))
            allocate(reo%recv_offset(sub_comm%size))

            allocate(reo%dest(box%lo(dir):box%hi(dir)))


            !--------------------------------------------------------------
            ! Obtain 'reverted' bounds of all MPI ranks

            allocate(rlos(0:layout%size(dir)-1))
            allocate(rhis(0:layout%size(dir)-1))
            do i = 0, layout%size(dir)-1
                call get_local_bounds(box%global_size(dir), i, layout%size(dir), lo, hi)
                rhis(i) = box%global_size(dir) - lo - 1
                rlos(i) = box%global_size(dir) - hi - 1
            enddo

            reo%send_recv_count = 0

            !--------------------------------------------------------------
            ! Determine destination rank and number of elements to send in x
            do i = box%lo(dir), box%hi(dir)
                reo%dest(i) = sub_comm%rank
                do j = 0, layout%size(dir)-1
                    if (i >= rlos(j) .and. i <= rhis(j)) then
                        call MPI_Cart_rank(sub_comm%comm, (/j/), rank, sub_comm%err)
                        reo%dest(i) = rank
                    endif
                enddo
                j = reo%dest(i)
                if (.not. sub_comm%rank == j) then
                    reo%send_recv_count(j+1) = reo%send_recv_count(j+1) + 1
                endif
            enddo

            !--------------------------------------------------------------
            ! Add y and z dimensions
            reo%recv_count = reo%send_recv_count
            reo%send_recv_count = reo%send_recv_count * box%size(d)

            deallocate(rlos)
            deallocate(rhis)

            allocate(reo%send_buffer(sum(reo%send_recv_count)))
            allocate(reo%recv_buffer(sum(reo%send_recv_count)))

            !--------------------------------------------------------------
            ! Determine offsets for send and receive buffers

            reo%send_offset(1) = 0
            do i = 2, sub_comm%size
                reo%send_offset(i) = sum(reo%send_recv_count(1:i-1))
                reo%recv_offset(sub_comm%size-i+1) = reo%send_offset(i)
            enddo

            reo%recv_offset(sub_comm%size) = 0
            do i = sub_comm%size-1, 1, -1
                reo%recv_offset(i) = sum(reo%send_recv_count(sub_comm%size:i+1:-1))
            enddo

            !--------------------------------------------------------------
            ! Allocate buffers for halo fill

            if (reo%dir == 1) then
                ! if in x direction
                allocate(reo%lo_buffer(box%size(2)))
                allocate(reo%hi_buffer(box%size(2)))
                allocate(reo%hi_halo_buffer(box%size(2)))
                allocate(reo%lo_halo_buffer(box%size(2)))
             else
                ! if in y direction
                allocate(reo%lo_buffer(box%size(1)))
                allocate(reo%hi_buffer(box%size(1)))
                allocate(reo%hi_halo_buffer(box%size(1)))
                allocate(reo%lo_halo_buffer(box%size(1)))
            endif

            call MPI_Cart_shift(sub_comm%comm, 0, 1, reo%lo_rank, reo%hi_rank, sub_comm%err)

        end subroutine initialise_reversing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_x(fs)
            double precision, intent(in) :: fs(box%lo(2):box%hi(2),   &
                                               box%hlo(1):box%hhi(1))
            integer                      :: ix, iy, j

            j = 1
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%lo(2), box%hi(2)
                    x_reo%send_buffer(j) = fs(iy, ix)
                    j = j + 1
                enddo
            enddo

        end subroutine copy_to_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_x(fs)
            double precision, intent(out) :: fs(box%lo(2):box%hi(2),   &
                                                box%hlo(1):box%hhi(1))
            integer                       :: ix, iy
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(2):box%hi(2))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(x_reo%recv_count)

            half_length = (box%size(1) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(1) + i
                k = box%hi(1) - n_recvs - i
                buf = fs(:, j)
                fs(:, j) = fs(:, k)
                fs(:, k) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(x_reo%recv_buffer)
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%hi(2), box%lo(2), -1
                    fs(iy, ix) = x_reo%recv_buffer(j)
                    j = j - 1
                enddo
            enddo

        end subroutine copy_from_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_y(fs)
            double precision, intent(in) :: fs(box%hlo(2):box%hhi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, j

            j = 1
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(iy)) then
                    cycle
                endif
                do ix = box%lo(1), box%hi(1)
                    y_reo%send_buffer(j) = fs(iy, ix)
                    j = j + 1
                enddo
            enddo

        end subroutine copy_to_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_y(fs)
            double precision, intent(out) :: fs(box%hlo(2):box%hhi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(1):box%hi(1))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(y_reo%recv_count)

            half_length = (box%size(2) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(2) + i
                k = box%hi(2) - n_recvs - i
                buf = fs(j, :)
                fs(j, :) = fs(k, :)
                fs(k, :) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(y_reo%recv_buffer)
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(iy)) then
                    cycle
                endif
                do ix = box%hi(1), box%lo(1), -1
                    fs(iy, ix) = y_reo%recv_buffer(j)
                    j = j - 1
                enddo
            enddo

        end subroutine copy_from_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reverse2d_x(fs, gs)
            double precision, intent(in)  :: fs(box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            double precision, intent(out) :: gs(box%lo(2):box%hi(2),   &
                                                box%hlo(1):box%hhi(1))
            integer                       :: send_size, recv_size

            gs(:, box%lo(1):box%hi(1)) = fs

            call copy_to_buffer_in_x(gs)

            send_size = size(x_reo%send_buffer)
            recv_size = size(x_reo%recv_buffer)

            call MPI_Alltoallv(x_reo%send_buffer(1:send_size),          &
                               x_reo%send_recv_count(1:x_comm%size),    &
                               x_reo%send_offset(1:x_comm%size),        &
                               MPI_DOUBLE_PRECISION,                    &
                               x_reo%recv_buffer(1:recv_size),          &
                               x_reo%send_recv_count(1:x_comm%size),    &
                               x_reo%recv_offset(1:x_comm%size),        &
                               MPI_DOUBLE_PRECISION,                    &
                               x_comm%comm,                             &
                               x_comm%err)

            call copy_from_buffer_in_x(gs)

            call halo_x_fill(gs)

        end subroutine reverse2d_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reverse2d_y(fs, gs)
            double precision, intent(in)  :: fs(box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision, intent(out) :: gs(box%hlo(2):box%hhi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: send_size, recv_size

            gs(box%lo(2):box%hi(2), :) = fs

            call copy_to_buffer_in_y(gs)

            send_size = size(y_reo%send_buffer)
            recv_size = size(y_reo%recv_buffer)

            call MPI_Alltoallv(y_reo%send_buffer(1:send_size),          &
                               y_reo%send_recv_count(1:y_comm%size),    &
                               y_reo%send_offset(1:y_comm%size),        &
                               MPI_DOUBLE_PRECISION,                    &
                               y_reo%recv_buffer(1:recv_size),          &
                               y_reo%send_recv_count(1:y_comm%size),    &
                               y_reo%recv_offset(1:y_comm%size),        &
                               MPI_DOUBLE_PRECISION,                    &
                               y_comm%comm,                             &
                               y_comm%err)

            call copy_from_buffer_in_y(gs)

            call halo_y_fill(gs)

        end subroutine reverse2d_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine halo_x_fill(gs)
            double precision, intent(inout) :: gs(box%lo(2):box%hi(2),   &
                                                  box%hlo(1):box%hhi(1))

            ! copy from interior to buffers
            x_reo%lo_buffer = gs(:, box%lo(1))
            x_reo%hi_buffer = gs(:, box%hi(1))

            call communicate_halo(x_reo, x_comm)

            ! copy from buffers to halo
            gs(:, box%hhi(1)) = x_reo%hi_halo_buffer
            gs(:, box%hlo(1)) = x_reo%lo_halo_buffer

        end subroutine halo_x_fill

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine halo_y_fill(gs)
            double precision, intent(inout) :: gs(box%hlo(2):box%hhi(2), &
                                                  box%lo(1):box%hi(1))

            ! copy from interior to buffers
            y_reo%lo_buffer = gs(box%lo(2), :)
            y_reo%hi_buffer = gs(box%hi(2), :)

            call communicate_halo(y_reo, y_comm)

            ! copy from buffers to halo
            gs(box%hhi(2), :) = y_reo%hi_halo_buffer
            gs(box%hlo(2), :) = y_reo%lo_halo_buffer

        end subroutine halo_y_fill

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine communicate_halo(reo, sub_comm)
            type(reorder_type), intent(inout) :: reo
            type(communicator), intent(inout) :: sub_comm
            type(MPI_Request)                 :: requests(2)
            type(MPI_Status)                  :: send_statuses(2)
            integer                           :: lb(1), ub(1)

            lb = lbound(reo%lo_buffer)
            ub = ubound(reo%lo_buffer)

            ! send west buffer to east halo
            call MPI_Isend(reo%lo_buffer(lb(1):ub(1)),  &
                           size(reo%lo_buffer),         &
                           MPI_DOUBLE_PRECISION,        &
                           reo%lo_rank,                 &
                           REVERSE_LO_TAG,              &
                           sub_comm%comm,               &
                           requests(1),                 &
                           sub_comm%err)

            call mpi_check_for_error(sub_comm, &
                "in MPI_Isend of mpi_reverse2d::communicate_halo.")

            lb = lbound(reo%hi_halo_buffer)
            ub = ubound(reo%hi_halo_buffer)

            ! receive west buffer to east halo (left to right)
            call MPI_Recv(reo%hi_halo_buffer(lb(1):ub(1)),  &
                          size(reo%hi_halo_buffer),         &
                          MPI_DOUBLE_PRECISION,             &
                          reo%hi_rank,                      &
                          REVERSE_LO_TAG,                   &
                          sub_comm%comm,                    &
                          MPI_STATUS_IGNORE,                &
                          sub_comm%err)

            call mpi_check_for_error(sub_comm, &
                "in MPI_Recv of mpi_reverse2d::communicate_halo.")

            lb = lbound(reo%hi_buffer)
            ub = ubound(reo%hi_buffer)

            ! send east buffer to west halo
            call MPI_Isend(reo%hi_buffer(lb(1):ub(1)),  &
                           size(reo%hi_buffer),         &
                           MPI_DOUBLE_PRECISION,        &
                           reo%hi_rank,                 &
                           REVERSE_HI_TAG,              &
                           sub_comm%comm,               &
                           requests(2),                 &
                           sub_comm%err)

            call mpi_check_for_error(sub_comm, &
                "in MPI_Isend of mpi_reverse2d::communicate_halo.")

            lb = lbound(reo%lo_halo_buffer)
            ub = ubound(reo%lo_halo_buffer)

            ! receive east buffer into west halo (right to left)
            call MPI_Recv(reo%lo_halo_buffer(lb(1):ub(1)),  &
                          size(reo%lo_halo_buffer),         &
                          MPI_DOUBLE_PRECISION,             &
                          reo%lo_rank,                      &
                          REVERSE_HI_TAG,                   &
                          sub_comm%comm,                    &
                          MPI_STATUS_IGNORE,                &
                          sub_comm%err)

            call mpi_check_for_error(sub_comm, &
                "in MPI_Recv of mpi_reverse2d::communicate_halo.")

            call MPI_Waitall(2,                 &
                            requests,           &
                            send_statuses,      &
                            sub_comm%err)

            call mpi_check_for_error(sub_comm, &
                "in MPI_Waitall of mpi_reverse2d::communicate_halo.")

        end subroutine communicate_halo

end module mpi_reverse2d
