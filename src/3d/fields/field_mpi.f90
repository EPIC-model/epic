module field_mpi
    use constants, only : zero
    use mpi_layout
    use mpi_communicator
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_message, mpi_check_for_error
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

            data(:, box%hlo(2):box%hhi(2),   box%hlo(1))              = zero  ! west halo
            data(:, box%hlo(2):box%hhi(2),   box%hhi(1)-1:box%hhi(1)) = zero  ! east halo
            data(:, box%hlo(2),              box%lo(1):box%hi(1))     = zero  ! south halo
            data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))     = zero  ! north halo

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
            double precision, dimension(:), pointer :: send_buf, recv_buf
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: statuses(8)
            integer                                 :: tag, n

            do n = 1, 8

                tag = NEIGHBOUR_TAG(n)

                call get_interior_buffer_ptr(tag, recv_buf)

                call MPI_Irecv(recv_buf,                &
                               size(recv_buf),          &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(tag)%rank,    &
                               tag,                     &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error("in MPI_Irecv of field_mpi::halo_to_interior_communication.")
            enddo

            do n = 1, 8

                call get_halo_buffer_ptr(n, send_buf)

                call MPI_Send(send_buf,                &
                              size(send_buf),          &
                              MPI_DOUBLE_PRECISION,    &
                              neighbours(n)%rank,      &
                              NEIGHBOUR_TAG(n),        &
                              comm%cart,               &
                              comm%err)

                call mpi_check_for_error("in MPI_Send of field_mpi::halo_to_interior_communication.")
            enddo

            call MPI_Waitall(8,         &
                             requests,  &
                             statuses,  &
                             comm%err)

            call mpi_check_for_error("in MPI_Waitall of field_mpi::halo_to_interior_communication.")

        end subroutine halo_to_interior_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine interior_to_halo_communication
            double precision, dimension(:), pointer :: send_buf, recv_buf
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: statuses(8)
            integer                                 :: tag, n

            do n = 1, 8

                tag = NEIGHBOUR_TAG(n)

                call get_halo_buffer_ptr(tag, recv_buf)

                call MPI_Irecv(recv_buf,                &
                               size(recv_buf),          &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(tag)%rank,    &
                               tag,                     &
                               comm%cart,               &
                               requests(n),             &
                               comm%err)

                call mpi_check_for_error("in MPI_Irecv of field_mpi::interior_to_halo_communication.")
            enddo


            do n = 1, 8
                call get_interior_buffer_ptr(n, send_buf)

                call MPI_Send(send_buf,                &
                              size(send_buf),          &
                              MPI_DOUBLE_PRECISION,    &
                              neighbours(n)%rank,      &
                              NEIGHBOUR_TAG(n),        &
                              comm%cart,               &
                              comm%err)

                call mpi_check_for_error("in MPI_Send of field_mpi::interior_to_halo_communication.")
            enddo

            call MPI_Waitall(8,         &
                             requests,  &
                             statuses,  &
                             comm%err)

            call mpi_check_for_error("in MPI_Waitall of field_mpi::interior_to_halo_communication.")

        end subroutine interior_to_halo_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine allocate_buffers
            integer :: zlen, ylen, xlen

            xlen = box%hi(1)-box%lo(1)+1
            ylen = box%hi(2)-box%lo(2)+1
            zlen = box%hhi(3)-box%hlo(3)+1

            allocate(west_buf(zlen, ylen, 2))
            allocate(east_halo_buf(zlen, ylen, 2))

            allocate(east_buf(zlen, ylen))
            allocate(west_halo_buf(zlen, ylen))

            allocate(south_buf(zlen, 2, xlen))
            allocate(north_halo_buf(zlen, 2, xlen))

            allocate(north_buf(zlen, xlen))
            allocate(south_halo_buf(zlen, xlen))

            allocate(northeast_buf(zlen))
            allocate(southwest_halo_buf(zlen))

            allocate(southeast_buf(zlen, 2))
            allocate(northwest_halo_buf(zlen, 2))

            allocate(southwest_buf(zlen, 2, 2))
            allocate(northeast_halo_buf(zlen, 2, 2))

            allocate(northwest_buf(zlen, 2))
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

            west_buf  = data(:, box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf  = data(:, box%lo(2):box%hi(2),   box%hi(1))
            south_buf = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf = data(:, box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf = data(:, box%hi(2),             box%hi(1))
            southeast_buf = data(:, box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf = data(:, box%hi(2),             box%lo(1):box%lo(1)+1)

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

            east_halo_buf = data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1))
            west_halo_buf = data(:, box%lo(2):box%hi(2), box%hlo(1))
            north_halo_buf = data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))
            south_halo_buf = data(:, box%hlo(2), box%lo(1):box%hi(1))
            southwest_halo_buf = data(:, box%hlo(2), box%hlo(1))
            northwest_halo_buf = data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1))
            northeast_halo_buf = data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1))
            southeast_halo_buf = data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1))

        end subroutine copy_from_halo_to_buffers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffers_to_halo(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add

            if (l_add) then

                data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) + east_halo_buf

                data(:, box%lo(2):box%hi(2), box%hlo(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hlo(1)) + west_halo_buf

                data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) + north_halo_buf

                data(:, box%hlo(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hlo(2), box%lo(1):box%hi(1)) + south_halo_buf

                data(:, box%hlo(2), box%hlo(1)) &
                    = data(:, box%hlo(2), box%hlo(1)) + southwest_halo_buf

                data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) + northwest_halo_buf

                data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    + northeast_halo_buf

                data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) + southeast_halo_buf
            else
                data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) = east_halo_buf
                data(:, box%lo(2):box%hi(2), box%hlo(1)) = west_halo_buf
                data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) = north_halo_buf
                data(:, box%hlo(2), box%lo(1):box%hi(1)) = south_halo_buf
                data(:, box%hlo(2), box%hlo(1)) = southwest_halo_buf
                data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) = northwest_halo_buf
                data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = northeast_halo_buf
                data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) = southeast_halo_buf
            endif

        end subroutine copy_from_buffers_to_halo

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine copy_from_buffers_to_interior(data, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_add
!             integer :: ylen, zlen, i, j, k

            if (l_add) then

!                     ylen = box%hi(2)-box%lo(2)+1
!                     zlen = box%hhi(3)-box%hlo(3)+1

!                     do i = 1, 2
!                         do j = 1, ylen
!                             do k = 1, zlen
!                                 print *, comm%rank, i, j, k, west_buf(k, j, i)
!                             enddo
!                         enddo
!                     enddo

!                     do i = box%lo(1), box%lo(1)+1
!                         do j = box%lo(2), box%hi(2)
!                             do k = box%hlo(3), box%hhi(3)
!                                 print *, comm%rank, i, j, k, data(k, j, i)
!                             enddo
! !                         enddo
!                     enddo

                data(:, box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(:, box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) + west_buf

                data(:, box%lo(2):box%hi(2), box%hi(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hi(1)) + east_buf

                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) &
                    = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) + south_buf

                data(:, box%hi(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hi(2), box%lo(1):box%hi(1)) + north_buf

                data(:, box%hi(2), box%hi(1)) &
                    = data(:, box%hi(2), box%hi(1)) + northeast_buf

                data(:, box%lo(2):box%lo(2)+1, box%hi(1)) &
                    = data(:, box%lo(2):box%lo(2)+1, box%hi(1)) + southeast_buf

                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) &
                    = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) + southwest_buf

                data(:, box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(:, box%hi(2), box%lo(1):box%lo(1)+1) + northwest_buf

            else
                data(:, box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1) = west_buf
                data(:, box%lo(2):box%hi(2),   box%hi(1))             = east_buf
                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))   = south_buf
                data(:, box%hi(2),             box%lo(1):box%hi(1))   = north_buf

                data(:, box%hi(2),             box%hi(1))             = northeast_buf
                data(:, box%lo(2):box%lo(2)+1, box%hi(1))             = southeast_buf
                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = southwest_buf
                data(:, box%hi(2),             box%lo(1):box%lo(1)+1) = northwest_buf
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

        ! See https://fortran-lang.org/en/learn/best_practices/multidim_arrays/
        ! about pointer with lower rank pointer to higher rank array
        ! (accessed 10 Feb 2023)
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

end module field_mpi
