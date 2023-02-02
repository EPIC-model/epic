! =============================================================================
!                       Test pencil module
!
!       This unit test checks the pencil transposition.
! =============================================================================
program test_mpi_pencil
    use unit_test
    use constants, only : zero, one, two, four
    use sta2dfft
    use pencil_fft
    use dimensions
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    implicit none

    double precision, allocatable :: values(:, :, :), vtrans(:, :, :)
    logical                       :: passed
    integer                       :: ix, iy !ix, iy, iz

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 8
    ny = 16
    nz = 32
    lower = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(vtrans(box%hlo(2):box%hhi(2), box%hlo(3):box%hhi(3), box%hlo(1):box%hhi(1)))

    values(:, :, :) = zero
    vtrans(:, :, :) = zero

!     print *, box%size
!     stop

!     do iz = box%lo(3), box%hi(3)
!         do ix = box%lo(1), box%hi(1)
            do iy = box%lo(2), box%hi(2)
                values(:, iy, :) = iy+1
            enddo
!             write(WRITE_VOR, '(1x,f13.6,6(1x,1p,e14.7))')
!             write(*, '(8(f1.0))')
!             print *, values(iz, 0:ny-1, ix)
!         enddo
!         print *, ""
!     enddo

    call initialise_pencil_fft(nx, ny, nz)


    call transpose_to_pencil(y_from_z_transposition, (/1, 2, 3/), dim_y_comm, FORWARD,        &
                             values(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                             fft_in_y_buffer)

    if (comm%rank == comm%master) then
!         do ix = box%lo(1), box%hi(1)
!             do iy = box%lo(2), box%hi(2)
                do iy = 1, ny
                    print *, fft_in_y_buffer(iy, 1, 1)
                enddo
!             enddo
!         enddo
    endif

    fft_in_y_buffer(:, :, :) = zero
    do ix = box%lo(1), box%hi(1)
        fft_in_y_buffer(:, ix+1-box%lo(1), :) = ix + 1
    enddo

    call transpose_to_pencil(x_from_y_transposition, (/2, 3, 1/), dim_x_comm, FORWARD, &
                             fft_in_y_buffer, fft_in_x_buffer)

    if (comm%rank == comm%master) then
        do ix = 1, nx
            print *, fft_in_x_buffer(ix, 1, 1)
        enddo
    endif

    fft_in_y_buffer(:, :, :) = zero
    fft_in_x_buffer(:, :, :) = zero
    print *, "rank", comm%rank, "size", size(fft_in_x_buffer, 3)
    do iy = 1, size(fft_in_x_buffer, 3)
        fft_in_x_buffer(:, :, iy) = iy
    enddo

    call transpose_to_pencil(y_from_x_transposition, (/3, 1, 2/), dim_x_comm, BACKWARD, &
                            fft_in_x_buffer, fft_in_y_buffer)

    call MPI_Barrier(comm%world)
    if (comm%rank == comm%master) then
        do iy = 1, ny
            print *, fft_in_y_buffer(iy, 1, 1)
        enddo
    endif

!     print *, "Z from X"
!
!     call transpose_to_pencil(z_from_x_transposition, (/3, 1, 2/), comm%cart, FORWARD, &
!                             fft_in_x_buffer,                                           &
!                             values(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
!
!     if (comm%rank == comm%master) then
!         do iy = box%lo(3), box%hi(3)
!             print *, values(iy, 1, 1)
!         enddo
!     endif

    call finalise_pencil_fft


    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test pencil transposition', passed)
    endif


end program test_mpi_pencil
