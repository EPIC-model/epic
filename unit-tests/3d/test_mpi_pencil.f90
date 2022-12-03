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
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    implicit none

    logical :: passed

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    nx = 16
    ny = 32
    nz = 64
    lower = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)



    call initialise_pencil_fft(nx, ny, nz)



    call finalise_pencil_fft


    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test pencil transposition', passed)
    endif


end program test_mpi_pencil
