module field_mpi
    use constants, only : zero
    use mpi_layout
    use mpi_environment
    use mpi_utils, only : mpi_exit_on_error, mpi_check_for_error, mpi_stop
    implicit none

        private

        double precision, contiguous, dimension(:, :, :, :), pointer :: west_buf,           &
                                                                        south_buf,          &
                                                                        southwest_buf,      &
                                                                        east_halo_buf,      &
                                                                        north_halo_buf,     &
                                                                        northeast_halo_buf


        double precision, contiguous, dimension(:, :, :), pointer :: east_buf,              &
                                                                     north_buf,             &
                                                                     southeast_buf,         &
                                                                     northwest_buf,         &
                                                                     west_halo_buf,         &
                                                                     south_halo_buf,        &
                                                                     northwest_halo_buf,    &
                                                                     southeast_halo_buf

        double precision, contiguous, dimension(:, :), pointer :: northeast_buf,        &
                                                                  southwest_halo_buf

        logical :: l_allocated = .false.

        integer :: n_comp = 0

        ! Convenient routines to fill halo cells for a single scalar field:
        interface field_halo_fill_scalar
            module procedure :: field_halo_fill_scalar_2d
            module procedure :: field_halo_fill_scalar_3d
        end interface field_halo_fill_scalar

        ! Convenient routines to fill halo cells for a single vector field:
        interface field_halo_fill_vector
            module procedure :: field_halo_fill_vector_2d
            module procedure :: field_halo_fill_vector_3d
        end interface field_halo_fill_vector

        ! Convenient routines to accumulate interior and halo cells
        ! for a single scalar field:
        interface field_halo_swap_scalar
            module procedure :: field_halo_swap_scalar_2d
            module procedure :: field_halo_swap_scalar_3d
        end interface field_halo_swap_scalar

        ! Convenient routines to accumulate interior and halo cells
        ! for a single vector field:
        interface field_halo_swap_vector
            module procedure :: field_halo_swap_vector_2d
            module procedure :: field_halo_swap_vector_3d
        end interface field_halo_swap_vector

        interface field_halo_reset
            module procedure :: field_halo_reset_scalar
            module procedure :: field_halo_reset_vector
        end interface field_halo_reset

        interface field_interior_accumulate_scalar
            module procedure :: field_interior_accumulate_scalar_2d
            module procedure :: field_interior_accumulate_scalar_3d
        end interface field_interior_accumulate_scalar

         interface field_interior_accumulate_vector
            module procedure :: field_interior_accumulate_vector_2d
            module procedure :: field_interior_accumulate_vector_3d
        end interface field_interior_accumulate_vector

        interface field_interior_to_buffer
            module procedure :: field_interior_to_buffer_2d
            module procedure :: field_interior_to_buffer_3d
        end interface field_interior_to_buffer

        interface field_buffer_to_halo
            module procedure :: field_buffer_to_halo_2d
            module procedure :: field_buffer_to_halo_3d
        end interface field_buffer_to_halo

        interface field_halo_to_buffer
            module procedure :: field_halo_to_buffer_2d
            module procedure :: field_halo_to_buffer_3d
        end interface field_halo_to_buffer

        interface field_buffer_to_interior
            module procedure :: field_buffer_to_interior_2d
            module procedure :: field_buffer_to_interior_3d
        end interface field_buffer_to_interior

        public :: field_halo_reset                  &
                , field_halo_fill_scalar            &
                , field_halo_fill_vector            &
                , field_interior_accumulate_scalar  &
                , field_interior_accumulate_vector  &
                , field_halo_swap_scalar            &
                , field_halo_swap_vector            &
                , field_mpi_alloc                   &
                , field_mpi_dealloc                 &
                , field_halo_to_buffer              &
                , field_buffer_to_halo              &
                , field_interior_to_buffer          &
                , field_buffer_to_interior          &
                , interior_to_halo_communication    &
                , halo_to_interior_communication    &
                , field_halo_to_buffer_integer      &
                , field_buffer_to_interior_integer  &
                , field_interior_to_buffer_integer  &
                , field_buffer_to_halo_integer

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_fill_scalar_3d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_alloc

            if (l_alloc) then
                call field_mpi_alloc(1, ndim=3)
            endif

            call field_interior_to_buffer_3d(data, 1)

            call interior_to_halo_communication

            call field_buffer_to_halo_3d(data, 1, .false.)

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_halo_fill_scalar_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_fill_scalar_2d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_alloc

            if (l_alloc) then
                call field_mpi_alloc(1, ndim=2)
            endif

            call field_interior_to_buffer_2d(data, 1)

            call interior_to_halo_communication

            call field_buffer_to_halo_2d(data, 1, .false.)

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_halo_fill_scalar_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_fill_vector_3d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(3):, &
                                                    box%hlo(2):, &
                                                    box%hlo(1):, &
                                                    :)
            logical,          intent(in)    :: l_alloc
            integer                         :: nc, ncomp

            ncomp = size(data, 4)

            if (l_alloc) then
                call field_mpi_alloc(ncomp, ndim=3)
            endif

            do nc = 1, ncomp
                call field_interior_to_buffer_3d(data(:, :, :, nc), nc)
            enddo

            call interior_to_halo_communication

            do nc = 1, ncomp
                call field_buffer_to_halo_3d(data(:, :, :, nc), nc, .false.)
            enddo

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_halo_fill_vector_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_fill_vector_2d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(2):, &
                                                    box%hlo(1):, &
                                                    :)
            logical,          intent(in)    :: l_alloc
            integer                         :: nc, ncomp

            ncomp = size(data, 3)

            if (l_alloc) then
                call field_mpi_alloc(ncomp, ndim=2)
            endif

            do nc = 1, ncomp
                call field_interior_to_buffer_2d(data(:, :, nc), nc)
            enddo

            call interior_to_halo_communication

            do nc = 1, ncomp
                call field_buffer_to_halo_2d(data(:, :, nc), nc, .false.)
            enddo

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_halo_fill_vector_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_reset_scalar(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))

            data(:, box%hlo(2):box%hhi(2),   box%hlo(1))              = zero  ! west halo
            data(:, box%hlo(2):box%hhi(2),   box%hhi(1)-1:box%hhi(1)) = zero  ! east halo
            data(:, box%hlo(2),              box%lo(1):box%hi(1))     = zero  ! south halo
            data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))     = zero  ! north halo

        end subroutine field_halo_reset_scalar

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_reset_vector(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1), &
                                                    3)

            data(:, box%hlo(2):box%hhi(2),   box%hlo(1), :)              = zero  ! west halo
            data(:, box%hlo(2):box%hhi(2),   box%hhi(1)-1:box%hhi(1), :) = zero  ! east halo
            data(:, box%hlo(2),              box%lo(1):box%hi(1), :)     = zero  ! south halo
            data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1), :)     = zero  ! north halo

        end subroutine field_halo_reset_vector

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_accumulate_scalar_3d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_alloc

            if (l_alloc) then
                call field_mpi_alloc(1, ndim=3)
            endif

            call field_halo_to_buffer_3d(data, 1)

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call field_buffer_to_interior_3d(data, 1, .true.)

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_interior_accumulate_scalar_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_accumulate_scalar_2d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            logical,          intent(in)    :: l_alloc

            if (l_alloc) then
                call field_mpi_alloc(1, ndim=2)
            endif

            call field_halo_to_buffer_2d(data, 1)

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call field_buffer_to_interior_2d(data, 1, .true.)

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_interior_accumulate_scalar_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_accumulate_vector_3d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(3):, &
                                                    box%hlo(2):, &
                                                    box%hlo(1):, &
                                                     :)
            logical,          intent(in)    :: l_alloc
            integer                         :: nc, ncomp

            ncomp = size(data, 4)

            if (l_alloc) then
                call field_mpi_alloc(ncomp, ndim=3)
            endif

            do nc = 1, ncomp
                call field_halo_to_buffer_3d(data(:, :, :, nc), nc)
            enddo

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            do nc = 1, ncomp
                call field_buffer_to_interior_3d(data(:, :, :, nc), nc, .true.)
            enddo

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_interior_accumulate_vector_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_accumulate_vector_2d(data, l_alloc)
            double precision, intent(inout) :: data(box%hlo(2):, &
                                                    box%hlo(1):, &
                                                     :)
            logical,          intent(in)    :: l_alloc
            integer                         :: nc, ncomp

            ncomp = size(data, 3)

            if (l_alloc) then
                call field_mpi_alloc(ncomp, ndim=2)
            endif

            do nc = 1, ncomp
                call field_halo_to_buffer_2d(data(:, :, nc), nc)
            enddo

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            do nc = 1, ncomp
                call field_buffer_to_interior_2d(data(:, :, nc), nc, .true.)
            enddo

            if (l_alloc) then
                call field_mpi_dealloc
            endif

        end subroutine field_interior_accumulate_vector_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_swap_scalar_3d(data)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_mpi_alloc(1, ndim=3)

            call field_interior_accumulate_scalar_3d(data, .false.)

            call field_halo_fill_scalar_3d(data, .false.)

            call field_mpi_dealloc

        end subroutine field_halo_swap_scalar_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_swap_scalar_2d(data)
            double precision, intent(inout) :: data(box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_mpi_alloc(1, ndim=2)

            call field_interior_accumulate_scalar_2d(data, .false.)

            call field_halo_fill_scalar_2d(data, .false.)

            call field_mpi_dealloc

        end subroutine field_halo_swap_scalar_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_swap_vector_3d(data)
            double precision, intent(inout) :: data(box%hlo(3):, &
                                                    box%hlo(2):, &
                                                    box%hlo(1):, &
                                                    :)
            integer                         :: ncomp
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.


            ncomp = size(data, 4)

            call field_mpi_alloc(ncomp, ndim=3)

            call field_interior_accumulate_vector_3d(data, .false.)

            call field_halo_fill_vector_3d(data, .false.)

            call field_mpi_dealloc

        end subroutine field_halo_swap_vector_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_swap_vector_2d(data)
            double precision, intent(inout) :: data(box%hlo(2):, &
                                                    box%hlo(1):, &
                                                    :)
            integer                         :: ncomp
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.


            ncomp = size(data, 3)

            call field_mpi_alloc(ncomp, ndim=2)

            call field_interior_accumulate_vector_2d(data, .false.)

            call field_halo_fill_vector_2d(data, .false.)

            call field_mpi_dealloc

        end subroutine field_halo_swap_vector_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine halo_to_interior_communication
            double precision, dimension(:), pointer :: send_buf, recv_buf
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: statuses(8)
            integer                                 :: tag, n, recv_size, send_size

            do n = 1, 8

                tag = RECV_NEIGHBOUR_TAG(n)

                call get_interior_buffer_ptr(tag, recv_buf)

                recv_size = size(recv_buf)

                call MPI_Irecv(recv_buf(1:recv_size),   &
                               recv_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               tag,                     &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Irecv of field_mpi::halo_to_interior_communication.")
            enddo

            do n = 1, 8

                call get_halo_buffer_ptr(n, send_buf)

                send_size = size(send_buf)

                call MPI_Send(send_buf(1:send_size),    &
                              send_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              SEND_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Send of field_mpi::halo_to_interior_communication.")
            enddo

            call MPI_Waitall(8,         &
                             requests,  &
                             statuses,  &
                             cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of field_mpi::halo_to_interior_communication.")

        end subroutine halo_to_interior_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine interior_to_halo_communication
            double precision, dimension(:), pointer :: send_buf, recv_buf
            type(MPI_Request)                       :: requests(8)
            type(MPI_Status)                        :: statuses(8)
            integer                                 :: tag, n, recv_size, send_size

            do n = 1, 8

                tag = RECV_NEIGHBOUR_TAG(n)

                call get_halo_buffer_ptr(tag, recv_buf)

                recv_size = size(recv_buf)

                call MPI_Irecv(recv_buf(1:recv_size),   &
                               recv_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               tag,                     &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Irecv of field_mpi::interior_to_halo_communication.")
            enddo


            do n = 1, 8
                call get_interior_buffer_ptr(n, send_buf)

                send_size = size(send_buf)

                call MPI_Send(send_buf(1:send_size),    &
                              send_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              SEND_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Send of field_mpi::interior_to_halo_communication.")
            enddo

            call MPI_Waitall(8,         &
                             requests,  &
                             statuses,  &
                             cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of field_mpi::interior_to_halo_communication.")

        end subroutine interior_to_halo_communication

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This routine is called inside fields::field_alloc
        ! @param[in] ncomp is the number of components (or number of scalar fields)
        ! @param[in] ndim is the number of dimensions. It can take either 2 or 3.
        subroutine field_mpi_alloc(ncomp, ndim)
            integer, intent(in) :: ncomp
            integer, intent(in) :: ndim
            integer             :: zlen, ylen, xlen

            if (l_allocated) then
                return
            endif
            l_allocated = .true.

            n_comp = ncomp

#ifndef NDEBUG
            if ((ndim /= 2) .and. (ndim /= 3)) then
                call mpi_stop("field_mpi::field_mpi_alloc can only allocate 2D or 3D buffers.")
            endif
#endif

            xlen = box%hi(1)-box%lo(1)+1
            ylen = box%hi(2)-box%lo(2)+1

            zlen = box%hhi(3)-box%hlo(3)+1
            if (ndim == 2) then
                zlen = 1
            endif

            allocate(west_buf(zlen, ylen, 2, n_comp))
            allocate(east_halo_buf(zlen, ylen, 2, n_comp))

            allocate(east_buf(zlen, ylen, n_comp))
            allocate(west_halo_buf(zlen, ylen, n_comp))

            allocate(south_buf(zlen, 2, xlen, n_comp))
            allocate(north_halo_buf(zlen, 2, xlen, n_comp))

            allocate(north_buf(zlen, xlen, n_comp))
            allocate(south_halo_buf(zlen, xlen, n_comp))

            allocate(northeast_buf(zlen, n_comp))
            allocate(southwest_halo_buf(zlen, n_comp))

            allocate(southeast_buf(zlen, 2, n_comp))
            allocate(northwest_halo_buf(zlen, 2, n_comp))

            allocate(southwest_buf(zlen, 2, 2, n_comp))
            allocate(northeast_halo_buf(zlen, 2, 2, n_comp))

            allocate(northwest_buf(zlen, 2, n_comp))
            allocate(southeast_halo_buf(zlen, 2, n_comp))

        end subroutine field_mpi_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This routine is called inside fields::field_dealloc
        subroutine field_mpi_dealloc

            if (.not. l_allocated) then
                return
            endif
            l_allocated = .false.
            n_comp = 0

            deallocate(west_buf)
            deallocate(east_halo_buf)

            deallocate(east_buf)
            deallocate(west_halo_buf)

            deallocate(south_buf)
            deallocate(north_halo_buf)

            deallocate(north_buf)
            deallocate(south_halo_buf)

            deallocate(northeast_buf)
            deallocate(southwest_halo_buf)

            deallocate(southeast_buf)
            deallocate(northwest_halo_buf)

            deallocate(southwest_buf)
            deallocate(northeast_halo_buf)

            deallocate(northwest_buf)
            deallocate(southeast_halo_buf)

        end subroutine field_mpi_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_to_buffer_integer(data, nc)
            integer, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                        box%hlo(2):box%hhi(2), &
                                        box%hlo(1):box%hhi(1))
            integer, intent(in) :: nc
            double precision    :: tmp(box%hlo(3):box%hhi(3), &
                                       box%hlo(2):box%hhi(2), &
                                       box%hlo(1):box%hhi(1))

            tmp = data

            call field_interior_to_buffer(tmp, nc)

        end subroutine field_interior_to_buffer_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_to_buffer_3d(data, nc)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))
            integer,          intent(in) :: nc

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_interior_to_buffer_3d: Component number outside bounds.")
            endif

            west_buf(:, :, :, nc)  = data(:, box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf(:, :, nc)     = data(:, box%lo(2):box%hi(2),   box%hi(1))
            south_buf(:, :, :, nc) = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf(:, :, nc)    = data(:, box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf(:, nc)       = data(:, box%hi(2),             box%hi(1))
            southeast_buf(:, :, nc)    = data(:, box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf(:, :, :, nc) = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf(:, :, nc)    = data(:, box%hi(2),             box%lo(1):box%lo(1)+1)

        end subroutine field_interior_to_buffer_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_interior_to_buffer_2d(data, nc)
            double precision, intent(in) :: data(box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))
            integer,          intent(in) :: nc

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_interior_to_buffer_2d: Component number outside bounds.")
            endif

            west_buf(1, :, :, nc)  = data(box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1)
            east_buf(1, :, nc)     = data(box%lo(2):box%hi(2),   box%hi(1))
            south_buf(1, :, :, nc) = data(box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))
            north_buf(1, :, nc)    = data(box%hi(2),             box%lo(1):box%hi(1))

            northeast_buf(1, nc)       = data(box%hi(2),             box%hi(1))
            southeast_buf(1, :, nc)    = data(box%lo(2):box%lo(2)+1, box%hi(1))
            southwest_buf(1, :, :, nc) = data(box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)
            northwest_buf(1, :, nc)    = data(box%hi(2),             box%lo(1):box%lo(1)+1)

        end subroutine field_interior_to_buffer_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_to_buffer_integer(data, nc)
            integer, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                        box%hlo(2):box%hhi(2), &
                                        box%hlo(1):box%hhi(1))
            integer, intent(in) :: nc
            double precision    :: tmp(box%hlo(3):box%hhi(3), &
                                       box%hlo(2):box%hhi(2), &
                                       box%hlo(1):box%hhi(1))

            !$omp parallel workshare
            tmp = data
            !$omp end parallel workshare

            call field_halo_to_buffer(tmp, nc)

        end subroutine field_halo_to_buffer_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine field_halo_to_buffer_3d(data, nc)
            double precision, intent(in) :: data(box%hlo(3):box%hhi(3), &
                                                 box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))
            integer,          intent(in) :: nc

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_halo_to_buffer_3d: Component number outside bounds.")
            endif

            east_halo_buf(:, :, :, nc)      = data(:, box%lo(2):box%hi(2),     box%hhi(1)-1:box%hhi(1))
            west_halo_buf(:, :, nc)         = data(:, box%lo(2):box%hi(2),     box%hlo(1))
            north_halo_buf(:, :, :, nc)     = data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))
            south_halo_buf(:, :, nc)        = data(:, box%hlo(2),              box%lo(1):box%hi(1))
            southwest_halo_buf(:, nc)       = data(:, box%hlo(2),              box%hlo(1))
            northwest_halo_buf(:, :, nc)    = data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1))
            northeast_halo_buf(:, :, :, nc) = data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1))
            southeast_halo_buf(:, :, nc)    = data(:, box%hlo(2),              box%hhi(1)-1:box%hhi(1))

        end subroutine field_halo_to_buffer_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_halo_to_buffer_2d(data, nc)
            double precision, intent(in) :: data(box%hlo(2):box%hhi(2), &
                                                 box%hlo(1):box%hhi(1))
            integer,          intent(in) :: nc

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_halo_to_buffer_2d: Component number outside bounds.")
            endif

            east_halo_buf(1, :, :, nc)      = data(box%lo(2):box%hi(2),     box%hhi(1)-1:box%hhi(1))
            west_halo_buf(1, :, nc)         = data(box%lo(2):box%hi(2),     box%hlo(1))
            north_halo_buf(1, :, :, nc)     = data(box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1))
            south_halo_buf(1, :, nc)        = data(box%hlo(2),              box%lo(1):box%hi(1))
            southwest_halo_buf(1, nc)       = data(box%hlo(2),              box%hlo(1))
            northwest_halo_buf(1, :, nc)    = data(box%hhi(2)-1:box%hhi(2), box%hlo(1))
            northeast_halo_buf(1, :, :, nc) = data(box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1))
            southeast_halo_buf(1, :, nc)    = data(box%hlo(2),              box%hhi(1)-1:box%hhi(1))

        end subroutine field_halo_to_buffer_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_halo_integer(data, nc, l_add)
            integer, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                           box%hlo(2):box%hhi(2), &
                                           box%hlo(1):box%hhi(1))
            integer, intent(in)    :: nc
            logical, intent(in)    :: l_add
            double precision       :: tmp(box%hlo(3):box%hhi(3), &
                                          box%hlo(2):box%hhi(2), &
                                          box%hlo(1):box%hhi(1))

            !$omp parallel workshare
            tmp = data
            !$omp end parallel workshare

            call field_buffer_to_halo(tmp, nc, l_add)

            !$omp parallel workshare
            data = int(tmp)
            !$omp end parallel workshare

        end subroutine field_buffer_to_halo_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_halo_3d(data, nc, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            integer,          intent(in)    :: nc
            logical,          intent(in)    :: l_add

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_buffer_to_halo_3d: Component number outside bounds.")
            endif

            if (l_add) then

                data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) + east_halo_buf(:, :, :, nc)

                data(:, box%lo(2):box%hi(2), box%hlo(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hlo(1)) + west_halo_buf(:, :, nc)

                data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) + north_halo_buf(:, :, :, nc)

                data(:, box%hlo(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hlo(2), box%lo(1):box%hi(1)) + south_halo_buf(:, :, nc)

                data(:, box%hlo(2), box%hlo(1)) &
                    = data(:, box%hlo(2), box%hlo(1)) + southwest_halo_buf(:, nc)

                data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) + northwest_halo_buf(:, :, nc)

                data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    + northeast_halo_buf(:, :, :, nc)

                data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) + southeast_halo_buf(:, :, nc)
            else
                data(:, box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) = east_halo_buf(:, :, :, nc)
                data(:, box%lo(2):box%hi(2), box%hlo(1)) = west_halo_buf(:, :, nc)
                data(:, box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) = north_halo_buf(:, :, :, nc)
                data(:, box%hlo(2), box%lo(1):box%hi(1)) = south_halo_buf(:, :, nc)
                data(:, box%hlo(2), box%hlo(1)) = southwest_halo_buf(:, nc)
                data(:, box%hhi(2)-1:box%hhi(2), box%hlo(1)) = northwest_halo_buf(:, :, nc)
                data(:, box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = northeast_halo_buf(:, :, :, nc)
                data(:, box%hlo(2), box%hhi(1)-1:box%hhi(1)) = southeast_halo_buf(:, :, nc)
            endif

        end subroutine field_buffer_to_halo_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_halo_2d(data, nc, l_add)
            double precision, intent(inout) :: data(box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            integer,          intent(in)    :: nc
            logical,          intent(in)    :: l_add

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_buffer_to_halo_2d: Component number outside bounds.")
            endif

            if (l_add) then

                data(box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) + east_halo_buf(1, :, :, nc)

                data(box%lo(2):box%hi(2), box%hlo(1)) &
                    = data(box%lo(2):box%hi(2), box%hlo(1)) + west_halo_buf(1, :, nc)

                data(box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) &
                    = data(box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) + north_halo_buf(1, :, :, nc)

                data(box%hlo(2), box%lo(1):box%hi(1)) &
                    = data(box%hlo(2), box%lo(1):box%hi(1)) + south_halo_buf(1, :, nc)

                data(box%hlo(2), box%hlo(1)) &
                    = data(box%hlo(2), box%hlo(1)) + southwest_halo_buf(1, nc)

                data(box%hhi(2)-1:box%hhi(2), box%hlo(1)) &
                    = data(box%hhi(2)-1:box%hhi(2), box%hlo(1)) + northwest_halo_buf(1, :, nc)

                data(box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) &
                    + northeast_halo_buf(1, :, :, nc)

                data(box%hlo(2), box%hhi(1)-1:box%hhi(1)) &
                    = data(box%hlo(2), box%hhi(1)-1:box%hhi(1)) + southeast_halo_buf(1, :, nc)
            else
                data(box%lo(2):box%hi(2), box%hhi(1)-1:box%hhi(1)) = east_halo_buf(1, :, :, nc)
                data(box%lo(2):box%hi(2), box%hlo(1)) = west_halo_buf(1, :, nc)
                data(box%hhi(2)-1:box%hhi(2), box%lo(1):box%hi(1)) = north_halo_buf(1, :, :, nc)
                data(box%hlo(2), box%lo(1):box%hi(1)) = south_halo_buf(1, :, nc)
                data(box%hlo(2), box%hlo(1)) = southwest_halo_buf(1, nc)
                data(box%hhi(2)-1:box%hhi(2), box%hlo(1)) = northwest_halo_buf(1, :, nc)
                data(box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = northeast_halo_buf(1, :, :, nc)
                data(box%hlo(2), box%hhi(1)-1:box%hhi(1)) = southeast_halo_buf(1, :, nc)
            endif

        end subroutine field_buffer_to_halo_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_interior_integer(data, nc, l_add)
            integer, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                           box%hlo(2):box%hhi(2), &
                                           box%hlo(1):box%hhi(1))
            integer, intent(in)    :: nc
            logical, intent(in)    :: l_add
            double precision       :: tmp(box%hlo(3):box%hhi(3), &
                                          box%hlo(2):box%hhi(2), &
                                          box%hlo(1):box%hhi(1))

            !$omp parallel workshare
            tmp = data
            !$omp end parallel workshare

            call field_buffer_to_interior(tmp, nc, l_add)

            !$omp parallel workshare
            data = int(tmp)
            !$omp end parallel workshare

        end subroutine field_buffer_to_interior_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_interior_3d(data, nc, l_add)
            double precision, intent(inout) :: data(box%hlo(3):box%hhi(3), &
                                                    box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            integer,          intent(in)    :: nc
            logical,          intent(in)    :: l_add

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_buffer_to_interior_3d: Component number outside bounds.")
            endif

            if (l_add) then

                data(:, box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(:, box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) + west_buf(:, :, :, nc)

                data(:, box%lo(2):box%hi(2), box%hi(1)) &
                    = data(:, box%lo(2):box%hi(2), box%hi(1)) + east_buf(:, :, nc)

                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) &
                    = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) + south_buf(:, :, :, nc)

                data(:, box%hi(2), box%lo(1):box%hi(1)) &
                    = data(:, box%hi(2), box%lo(1):box%hi(1)) + north_buf(:, :, nc)

                data(:, box%hi(2), box%hi(1)) &
                    = data(:, box%hi(2), box%hi(1)) + northeast_buf(:, nc)

                data(:, box%lo(2):box%lo(2)+1, box%hi(1)) &
                    = data(:, box%lo(2):box%lo(2)+1, box%hi(1)) + southeast_buf(:, :, nc)

                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) &
                    = data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) + southwest_buf(:, :, :, nc)

                data(:, box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(:, box%hi(2), box%lo(1):box%lo(1)+1) + northwest_buf(:, :, nc)

            else
                data(:, box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1) = west_buf(:, :, :, nc)
                data(:, box%lo(2):box%hi(2),   box%hi(1))             = east_buf(:, :, nc)
                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))   = south_buf(:, :, :, nc)
                data(:, box%hi(2),             box%lo(1):box%hi(1))   = north_buf(:, :, nc)

                data(:, box%hi(2),             box%hi(1))             = northeast_buf(:, nc)
                data(:, box%lo(2):box%lo(2)+1, box%hi(1))             = southeast_buf(:, :, nc)
                data(:, box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = southwest_buf(:, :, :, nc)
                data(:, box%hi(2),             box%lo(1):box%lo(1)+1) = northwest_buf(:, :, nc)
            endif

        end subroutine field_buffer_to_interior_3d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_buffer_to_interior_2d(data, nc, l_add)
            double precision, intent(inout) :: data(box%hlo(2):box%hhi(2), &
                                                    box%hlo(1):box%hhi(1))
            integer,          intent(in)    :: nc
            logical,          intent(in)    :: l_add

            if ((nc < 1) .or. (nc > n_comp)) then
                call mpi_exit_on_error(&
                    "field_mpi::field_buffer_to_interior_2d: Component number outside bounds.")
            endif

            if (l_add) then

                data(box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%lo(2):box%hi(2), box%lo(1):box%lo(1)+1) + west_buf(1, :, :, nc)

                data(box%lo(2):box%hi(2), box%hi(1)) &
                    = data(box%lo(2):box%hi(2), box%hi(1)) + east_buf(1, :, nc)

                data(box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) &
                    = data(box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1)) + south_buf(1, :, :, nc)

                data(box%hi(2), box%lo(1):box%hi(1)) &
                    = data(box%hi(2), box%lo(1):box%hi(1)) + north_buf(1, :, nc)

                data(box%hi(2), box%hi(1)) &
                    = data(box%hi(2), box%hi(1)) + northeast_buf(1, nc)

                data(box%lo(2):box%lo(2)+1, box%hi(1)) &
                    = data(box%lo(2):box%lo(2)+1, box%hi(1)) + southeast_buf(1, :, nc)

                data(box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) &
                    = data(box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) + southwest_buf(1, :, :, nc)

                data(box%hi(2), box%lo(1):box%lo(1)+1) &
                    = data(box%hi(2), box%lo(1):box%lo(1)+1) + northwest_buf(1, :, nc)

            else
                data(box%lo(2):box%hi(2),   box%lo(1):box%lo(1)+1) = west_buf(1, :, :, nc)
                data(box%lo(2):box%hi(2),   box%hi(1))             = east_buf(1, :, nc)
                data(box%lo(2):box%lo(2)+1, box%lo(1):box%hi(1))   = south_buf(1, :, :, nc)
                data(box%hi(2),             box%lo(1):box%hi(1))   = north_buf(1, :, nc)

                data(box%hi(2),             box%hi(1))             = northeast_buf(1, nc)
                data(box%lo(2):box%lo(2)+1, box%hi(1))             = southeast_buf(1, :, nc)
                data(box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = southwest_buf(1, :, :, nc)
                data(box%hi(2),             box%lo(1):box%lo(1)+1) = northwest_buf(1, :, nc)
            endif

        end subroutine field_buffer_to_interior_2d

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
