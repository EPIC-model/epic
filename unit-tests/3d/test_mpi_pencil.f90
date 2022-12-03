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

!     double precision              :: error = zero


    nx = 16
    ny = 32
    nz = 64
    lower = (/zero, zero, zero/)
    extent = (/one, two, four/)

    call update_parameters

    call initialise_pencil_fft(nx, ny, nz)



    call finalise_pencil_fft


end program test_mpi_pencil
