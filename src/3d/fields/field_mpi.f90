module field_mpi
    use constants, only : zero
    use mpi_layout
    use mpi_communicator
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_message
    implicit none

        double precision, contiguous, dimension(:, :, :), pointer :: west_buf,           &
                                                                                south_buf,          &
                                                                                southwest_buf,      &
                                                                                east_halo_buf,      &
                                                                                north_halo_buf,     &
                                                                                northeast_halo_buf


        double precision, contiguous, dimension(:, :), pointer :: east_buf,             &
                                                                              north_buf,            &
                                                                              southeast_buf,        &
                                                                              northwest_buf,        &
                                                                              west_halo_buf,        &
                                                                              south_halo_buf,       &
                                                                              northwest_halo_buf,   &
                                                                              southeast_halo_buf

        double precision, contiguous, dimension(:), pointer :: northeast_buf,       &
                                                                           southwest_halo_buf

        logical :: l_allocated = .false.

        public :: field_halo_reset,          &
                  field_halo_fill,           &
                  field_interior_accumulate, &
                  field_halo_swap

    contains

        subroutine field_halo_reset(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            data(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2),   box%hlo(1))              = zero  ! west halo
            data(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2),   box%hhi(1)-1:box%hhi(1)) = zero  ! east halo
            data(box%hlo(3):box%hhi(3), box%hlo(2),              box%lo(1):box%hi(1))     = zero  ! south halo
            data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))     = zero  ! north halo

        end subroutine field_halo_reset

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_fill(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            call copy_from_interior_to_buffers(data)

            call interior_to_halo_communication

            call copy_from_buffers_to_halo(data, .false.)

        end subroutine field_halo_fill

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_accumulate(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            call copy_from_halo_to_buffers(data)

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call copy_from_buffers_to_interior(data, .true.)

        end subroutine field_interior_accumulate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_swap(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_interior_accumulate(data)

            call field_halo_fill(data)

        end subroutine field_halo_swap

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine halo_to_interior_communication
            double precision, dimension(:), pointer :: sendbuf, recvbuf
            integer                                 :: send_size, recv_size
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: tag, source, nb, n

            do n = 1, 8
                call get_halo_buffer_ptr(n, sendbuf)

                send_size = size(sendbuf)

                call MPI_Isend(sendbuf,                 &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                nb = map_tag_to_neighbour(tag)

                call get_interior_buffer_ptr(nb, recvbuf)

                call MPI_Recv(recvbuf,                  &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

        end subroutine halo_to_interior_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine interior_to_halo_communication
            double precision, dimension(:), pointer :: sendbuf, recvbuf
            integer                                 :: send_size, recv_size
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: recv_status, send_statuses(8)
            integer                                 :: tag, source, nb, n

            do n = 1, 8
                call get_interior_buffer_ptr(n, sendbuf)

                send_size = size(sendbuf)

                call MPI_Isend(sendbuf,                 &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               NEIGHBOUR_TAG(n),        &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)
            enddo

            do n = 1, 8

                ! check for incoming messages
                call mpi_check_for_message(tag, recv_size, source)

                nb = map_tag_to_neighbour(tag)

                call get_halo_buffer_ptr(nb, recvbuf)

                call MPI_Recv(recvbuf,                  &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              source,                   &
                              tag,                      &
                              comm%cart,                &
                              recv_status,              &
                              comm%err)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            comm%err)

        end subroutine interior_to_halo_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_buffers
            integer :: zlen, ylen, xlen

            xlen = box%hi(1)-box%lo(1)+1
            ylen = box%hi(2)-box%lo(2)+1
            zlen = box%hhi(3)-box%hlo(3)+1

            ! in the context of the sender process
            allocate(west_buf(zlen, ylen, 2))
            allocate(east_buf(zlen, ylen))
            allocate(south_buf(zlen, 2, xlen))
            allocate(north_buf(zlen, xlen))

            allocate(northeast_buf(zlen))
            allocate(southeast_buf(zlen, 2))
            allocate(southwest_buf(zlen, 2, 2))
            allocate(northwest_buf(zlen, 2))

            allocate(east_halo_buf(zlen, ylen, 2))
            allocate(west_halo_buf(zlen, ylen))
            allocate(north_halo_buf(zlen, 2, xlen))
            allocate(south_halo_buf(zlen, xlen))

            allocate(southwest_halo_buf(zlen))
            allocate(northwest_halo_buf(zlen, 2))
            allocate(northeast_halo_buf(zlen, 2, 2))
            allocate(southeast_halo_buf(zlen, 2))

        end subroutine allocate_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_interior_to_buffers(data)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))

            if (.not. l_allocated) then
                l_allocated = .true.
                call allocate_buffers
            endif

            west_buf  = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf  = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2),   box%hi(1))
            south_buf = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf = data(box%hlo(3):box%hhi(3), box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf = data(box%hlo(3):box%hhi(3), box%hi(2),             box%hi(1))
            southeast_buf = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf = data(box%hlo(3):box%hhi(3), box%hi(2),             box%lo(1):box%lo(1)+1)

        end subroutine copy_from_interior_to_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_halo_to_buffers(data)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))

            if (.not. l_allocated) then
                l_allocated = .true.
                call allocate_buffers
            endif

            east_halo_buf = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1))
            west_halo_buf = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hlo(1))
            north_halo_buf = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))
            south_halo_buf = data(box%hlo(3):box%hhi(3), box%hlo(2), box%lo(1):box%hi(1))
            southwest_halo_buf = data(box%hlo(3):box%hhi(3), box%hlo(2), box%hlo(1))
            northwest_halo_buf = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1))
            northeast_halo_buf = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1))
            southeast_halo_buf = data(box%hlo(3):box%hhi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1))

        end subroutine copy_from_halo_to_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffers_to_halo(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add

            if (l_add) then

                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) + east_halo_buf

                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hlo(1)) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hlo(1)) + west_halo_buf

                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) + north_halo_buf

                data(box%hlo(3):box%hhi(3), box%hlo(2), box%lo(1):box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hlo(2), box%lo(1):box%hi(1)) + south_halo_buf

                data(box%hlo(3):box%hhi(3), box%hlo(2), box%hlo(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hlo(2), box%hlo(1)) + southwest_halo_buf

                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) + northwest_halo_buf

                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    + northeast_halo_buf

                data(box%hlo(3):box%hhi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) + southeast_halo_buf
            else
                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) = east_halo_buf
                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hlo(1)) = west_halo_buf
                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) = north_halo_buf
                data(box%hlo(3):box%hhi(3), box%hlo(2), box%lo(1):box%hi(1)) = south_halo_buf
                data(box%hlo(3):box%hhi(3), box%hlo(2), box%hlo(1)) = southwest_halo_buf
                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) = northwest_halo_buf
                data(box%hlo(3):box%hhi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = northeast_halo_buf
                data(box%hlo(3):box%hhi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) = southeast_halo_buf
            endif

        end subroutine copy_from_buffers_to_halo

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffers_to_interior(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add

            if (l_add) then
                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) + west_buf

                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2), box%hi(1)) + east_buf

                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) + south_buf

                data(box%hlo(3):box%hhi(3), box%hi(2), box%lo(1):box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hi(2), box%lo(1):box%hi(1)) + north_buf

                data(box%hlo(3):box%hhi(3), box%hi(2), box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%hi(2), box%hi(1)) + northeast_buf

                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%hi(1)) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%hi(1)) + southeast_buf

                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) &
                    = data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) + southwest_buf

                data(box%hlo(3):box%hhi(3), box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%hlo(3):box%hhi(3), box%hi(2), box%lo(1):box%lo(1)+1) + northwest_buf

            else
                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1) = west_buf
                data(box%hlo(3):box%hhi(3), box%lo(2):box%hi(2),   box%hi(1))             = east_buf
                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))   = south_buf
                data(box%hlo(3):box%hhi(3), box%hi(2),             box%lo(1):box%hi(1))   = north_buf

                data(box%hlo(3):box%hhi(3), box%hi(2),             box%hi(1))             = northeast_buf
                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%hi(1))             = southeast_buf
                data(box%hlo(3):box%hhi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = southwest_buf
                data(box%hlo(3):box%hhi(3), box%hi(2),             box%lo(1):box%lo(1)+1) = northwest_buf
            endif

        end subroutine copy_from_buffers_to_interior

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_interior_buffer_ptr(dir, buf_ptr)
            integer,                                 intent(in)  :: dir
            double precision, dimension(:), pointer, intent(out) :: buf_ptr

            select case (dir)
                case (MPI_NORTH)
                    buf_ptr(1:size(north_buf)) => north_buf
                case (MPI_SOUTH)
                    buf_ptr(1:size(south_buf)) => south_buf
                case (MPI_WEST)
                    buf_ptr(1:size(west_buf)) => west_buf
                case (MPI_EAST)
                    buf_ptr(1:size(east_buf)) => east_buf
                case (MPI_NORTHWEST)
                    buf_ptr(1:size(northwest_buf)) => northwest_buf
                case (MPI_NORTHEAST)
                    buf_ptr(1:size(northeast_buf)) => northeast_buf
                case (MPI_SOUTHWEST)
                    buf_ptr(1:size(southwest_buf)) => southwest_buf
                case (MPI_SOUTHEAST)
                    buf_ptr(1:size(southeast_buf)) => southeast_buf
                case default
                    call mpi_exit_on_error("get_interior_buffer_ptr: No valid direction.")
            end select

        end subroutine get_interior_buffer_ptr

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine get_halo_buffer_ptr(dir, buf_ptr)
            integer,                                 intent(in)  :: dir
            double precision, dimension(:), pointer, intent(out) :: buf_ptr

            select case (dir)
                case (MPI_NORTH)
                    buf_ptr(1:size(north_halo_buf)) => north_halo_buf
                case (MPI_SOUTH)
                    buf_ptr(1:size(south_halo_buf)) => south_halo_buf
                case (MPI_WEST)
                    buf_ptr(1:size(west_halo_buf)) => west_halo_buf
                case (MPI_EAST)
                    buf_ptr(1:size(east_halo_buf)) => east_halo_buf
                case (MPI_NORTHWEST)
                    buf_ptr(1:size(northwest_halo_buf)) => northwest_halo_buf
                case (MPI_NORTHEAST)
                    buf_ptr(1:size(northeast_halo_buf)) => northeast_halo_buf
                case (MPI_SOUTHWEST)
                    buf_ptr(1:size(southwest_halo_buf)) => southwest_halo_buf
                case (MPI_SOUTHEAST)
                    buf_ptr(1:size(southeast_halo_buf)) => southeast_halo_buf
                case default
                    call mpi_exit_on_error("get_halo_buffer_ptr: No valid direction.")
            end select

        end subroutine get_halo_buffer_ptr

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Example: If the north neighbour sends, then we must return south.
        function map_tag_to_neighbour(tag) result(nb)
            integer, intent(in) :: tag
            integer             :: nb

            select case (tag)
                case (MPI_NORTH)
                    nb = MPI_SOUTH
                case (MPI_SOUTH)
                    nb = MPI_NORTH
                case (MPI_WEST)
                    nb = MPI_EAST
                case (MPI_EAST)
                    nb = MPI_WEST
                case (MPI_NORTHWEST)
                    nb = MPI_SOUTHEAST
                case (MPI_NORTHEAST)
                    nb = MPI_SOUTHWEST
                case (MPI_SOUTHWEST)
                    nb = MPI_NORTHEAST
                case (MPI_SOUTHEAST)
                    nb = MPI_NORTHWEST
                case default
                    call mpi_exit_on_error("map_tag_to_neighbour: No valid neighbour tag.")
            end select

        end function map_tag_to_neighbour

end module field_mpi
