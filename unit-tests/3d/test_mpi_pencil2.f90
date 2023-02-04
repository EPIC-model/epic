! =============================================================================
!                       Test pencil module
!
!       This unit test checks the pencil transposition.
! =============================================================================
program test_mpi_pencil2
    use unit_test
    use constants, only : zero, one, two, four
    use sta2dfft
    use fft_pencil
    use dimensions
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    implicit none

    double precision, allocatable :: values(:, :, :), vcheck(:, :, :)
    logical                       :: passed
    integer                       :: ix, iy, iz

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
    allocate(vcheck(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    values(:, :, :) = zero

    do iz = box%lo(3), box%hi(3)
        do iy = box%lo(2), box%hi(2)
            do ix = box%lo(1), box%hi(1)
                values(iz, iy, ix) = 1 + ix + nx*iy + nx*ny*iz
            enddo
        enddo
    enddo

    vcheck = values

    call initialise_pencil_fft(nx, ny, nz)


    call transpose_to_pencil(y_from_z_transposition, (/1, 2, 3/), dim_y_comm, FORWARD,        &
                             values(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                             fft_in_y_buffer)

    call transpose_to_pencil(x_from_y_transposition, (/2, 3, 1/), dim_x_comm, FORWARD, &
                             fft_in_y_buffer, fft_in_x_buffer)

    call transpose_to_pencil(y_from_x_transposition, (/3, 1, 2/), dim_x_comm, BACKWARD, &
                            fft_in_x_buffer, fft_in_y_buffer)

    values(:, :, :) = zero
    call transpose_to_pencil(z_from_y_transposition, (/2, 3, 1/), dim_y_comm, BACKWARD, &
                            fft_in_y_buffer,                                            &
                            values(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call MPI_Barrier(comm%world)

    print *, "rank", comm%rank, "check:", (0.0 == sum(dabs(values - vcheck)))


    call finalise_pencil_fft


    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI pencil 2 transposition', passed)
    endif


end program test_mpi_pencil2
