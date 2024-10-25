! =============================================================================
!                         Test MPI field diagnostics
! =============================================================================
program test_mpi_field_diagnostics
    use constants, only : zero, one, f12
    use unit_test
    use mpi_environment
    use mpi_layout
    use fields
    use field_diagnostics
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, vcell, vcelli, ncell
    use mpi_timer
    implicit none

    logical :: passed = .true.

    call mpi_env_initialise

    passed = (world%err == 0)

    call register_timer('field stats', field_stats_timer)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call field_alloc

    volg = vcell + one
    nparg = 1
    nsparg = 1
    velog = f12

    call calculate_field_diagnostics

    if (world%rank == world%root) then
        passed = (passed .and. (abs(field_stats(IDX_RMS_V) - vcelli) == zero))
        passed = (passed .and. (abs(field_stats(IDX_ABSERR_V) - vcelli) == zero))
        passed = (passed .and. (field_stats(IDX_MAX_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_MIN_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_AVG_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_AVG_NSPAR) == one))
        passed = (passed .and. abs(field_stats(IDX_KEG) - 0.375d0 * dble(ncell) * (vcell + one)) < 1.0e-14)
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI field diagnostics', passed)
    endif

end program test_mpi_field_diagnostics
