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
    integer                       :: ix, iy, iz

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    nx = 16
    ny = 32
    nz = 64
    lower = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(vtrans(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1), box%hlo(3):box%hhi(3)))

    values(:, :, :) = zero
    vtrans(:, :, :) = zero

    do ix = box%lo(1), box%hi(1)
        do iy = box%lo(2), box%hi(2)
            do iz = box%lo(3), box%hi(3)
                values(iz, iy, ix) = iy+1
            enddo
        enddo
    enddo


    call initialise_pencil_fft(nx, ny, nz)


    call transpose_to_pencil(y_from_z_transposition, (/I_Z, I_Y, I_X/), dim_y_comm, FORWARD, values, vtrans)

    if (mpi_rank == mpi_master) then
!         do ix = box%lo(1), box%hi(1)
!             do iy = box%lo(2), box%hi(2)
                do iy = box%lo(2), box%hi(2)
                    print *, vtrans(iy, box%lo(1), box%lo(3))
                enddo
!             enddo
!         enddo
    endif

    call finalise_pencil_fft


    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test pencil transposition', passed)
    endif


end program test_mpi_pencil
