module mpi_reverse
    use mpi_communicator
    use mpi_layout
    use mpi_utils, only : mpi_exit_on_error
    use mpi_halo, only : halo_fill
    implicit none

    type :: reorder_type
        integer                       :: dir                            ! reorder index
        integer, allocatable          :: dest(:)                        ! destination rank
        integer, allocatable          :: send_recv_count(:)
        integer, allocatable          :: recv_count(:)                  ! number of receives of
                                                                        ! the reversing dimension
        integer, allocatable          :: send_offset(:), recv_offset(:)
        double precision, allocatable :: send_buffer(:), recv_buffer(:)
    end type reorder_type

    logical :: l_initialised_x = .false.
    logical :: l_initialised_y = .false.

    private :: copy_from_buffer_in_x,   &
               copy_to_buffer_in_x,     &
               copy_from_buffer_in_y,   &
               copy_to_buffer_in_y,     &
               l_initialised_x,         &
               l_initialised_y,         &
               setup_reversing

    type(sub_communicator) :: x_comm
    type(sub_communicator) :: y_comm

    type(reorder_type) :: x_reo
    type(reorder_type) :: y_reo

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine setup_reversing(reo, sub_comm, dir)
            type(reorder_type),     intent(inout) :: reo
            type(sub_communicator), intent(inout) :: sub_comm
            integer,                intent(in)    :: dir
            integer, allocatable                  :: rlos(:), rhis(:)
            integer                               :: i, j, lo, hi, rank, d
            logical                               :: remain_dims(2)

            reo%dir = dir

            ! if dir = 1 --> d = 2, i.e. if dir = x --> d = y
            ! if dir = 2 --> d = 1, i.e. if dir = y --> d = y
            d = mod(dir, 2) + 1

            remain_dims(dir) = .true.
            remain_dims(d) = .false.

            call MPI_Cart_sub(comm%cart, remain_dims, sub_comm%comm, comm%err)

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
                call get_local_bounds(box%global_ncells(dir), i, layout%size(dir), lo, hi)
                rhis(i) = box%global_ncells(dir) - lo - 1
                rlos(i) = box%global_ncells(dir) - hi - 1
            enddo

            reo%send_recv_count = 0

            !--------------------------------------------------------------
            ! Determine destination rank and number of elements to send in x
            do i = box%lo(dir), box%hi(dir)
                reo%dest(i) = sub_comm%rank
                do j = 0, layout%size(dir)-1
                    if (i >= rlos(j) .and. i <= rhis(j)) then
                        call MPI_Cart_rank(sub_comm%comm, (/j/), rank, comm%err)
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
            reo%send_recv_count = reo%send_recv_count * box%size(d) * box%size(3)

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

        end subroutine setup_reversing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_x(fs)
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz, j

            j = 1
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%lo(2), box%hi(2)
                    do iz = box%lo(3), box%hi(3)
                        x_reo%send_buffer(j) = fs(iz, iy, ix)
                        j = j + 1
                    enddo
                enddo
            enddo

        end subroutine copy_to_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_x(fs)
            double precision, intent(out) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy, iz
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(3):box%hi(3), box%lo(2):box%hi(2))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(x_reo%recv_count)

            half_length = (box%size(1) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(1) + i
                k = box%hi(1) - n_recvs - i
                buf = fs(:, :, j)
                fs(:, :, j) = fs(:, :, k)
                fs(:, :, k) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(x_reo%recv_buffer)
            do ix = box%hi(1), box%lo(1), -1
                if (x_comm%rank == x_reo%dest(ix)) then
                    cycle
                endif
                do iy = box%hi(2), box%lo(2), -1
                    do iz = box%hi(3), box%lo(3), -1
                        fs(iz, iy, ix) = x_reo%recv_buffer(j)
                        j = j - 1
                    enddo
                enddo
            enddo

        end subroutine copy_from_buffer_in_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_to_buffer_in_y(fs)
            double precision, intent(in) :: fs(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1))
            integer                      :: ix, iy, iz, j

            j = 1
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(iy)) then
                    cycle
                endif
                do ix = box%lo(1), box%hi(1)
                    do iz = box%lo(3), box%hi(3)
                        y_reo%send_buffer(j) = fs(iz, iy, ix)
                        j = j + 1
                    enddo
                enddo
            enddo

        end subroutine copy_to_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffer_in_y(fs)
            double precision, intent(out) :: fs(box%lo(3):box%hi(3), &
                                                box%lo(2):box%hi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ix, iy, iz
            integer                       :: i, j, k, half_length, n_recvs
            double precision              :: buf(box%lo(3):box%hi(3), box%lo(1):box%hi(1))


            !--------------------------------------------------------------
            ! Revert array locally

            n_recvs = sum(y_reo%recv_count)

            half_length = (box%size(2) - n_recvs) / 2
            do i = 0, half_length-1
                j = box%lo(2) + i
                k = box%hi(2) - n_recvs - i
                buf = fs(:, j, :)
                fs(:, j, :) = fs(:, k, :)
                fs(:, k, :) = buf
            enddo

            !--------------------------------------------------------------
            ! Copy from buffer to actual array

            j = size(y_reo%recv_buffer)
            do iy = box%hi(2), box%lo(2), -1
                if (y_comm%rank == y_reo%dest(iy)) then
                    cycle
                endif
                do ix = box%hi(1), box%lo(1), -1
                    do iz = box%hi(3), box%lo(3), -1
                        fs(iz, iy, ix) = y_reo%recv_buffer(j)
                        j = j - 1
                    enddo
                enddo
            enddo

        end subroutine copy_from_buffer_in_y

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reverse_x(fs, gs)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))

            if (.not. l_initialised_x) then
                call setup_reversing(x_reo, x_comm, 1)
                l_initialised_x = .true.
            endif

            gs = fs

            call copy_to_buffer_in_x(fs(box%lo(3):box%hi(3), &
                                        box%lo(2):box%hi(2), &
                                        box%lo(1):box%hi(1)))

            call MPI_alltoallv(x_reo%send_buffer,       &
                               x_reo%send_recv_count,   &
                               x_reo%send_offset,       &
                               MPI_DOUBLE_PRECISION,    &
                               x_reo%recv_buffer,       &
                               x_reo%send_recv_count,   &
                               x_reo%recv_offset,       &
                               MPI_DOUBLE_PRECISION,    &
                               x_comm%comm,             &
                               x_comm%err)

            call copy_from_buffer_in_x(gs(box%lo(3):box%hi(3), &
                                          box%lo(2):box%hi(2), &
                                          box%lo(1):box%hi(1)))

            call halo_fill(gs)

        end subroutine reverse_x

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reverse_y(fs, gs)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))

            if (.not. l_initialised_y) then
                call setup_reversing(y_reo, y_comm, 2)
                l_initialised_y = .true.
            endif

            gs = fs

            call copy_to_buffer_in_y(fs(box%lo(3):box%hi(3), &
                                        box%lo(2):box%hi(2), &
                                        box%lo(1):box%hi(1)))

            call MPI_alltoallv(y_reo%send_buffer,       &
                               y_reo%send_recv_count,   &
                               y_reo%send_offset,       &
                               MPI_DOUBLE_PRECISION,    &
                               y_reo%recv_buffer,       &
                               y_reo%send_recv_count,   &
                               y_reo%recv_offset,       &
                               MPI_DOUBLE_PRECISION,    &
                               y_comm%comm,             &
                               y_comm%err)

            call copy_from_buffer_in_y(gs(box%lo(3):box%hi(3), &
                                          box%lo(2):box%hi(2), &
                                          box%lo(1):box%hi(1)))

            call halo_fill(gs)

        end subroutine reverse_y

end module mpi_reverse
