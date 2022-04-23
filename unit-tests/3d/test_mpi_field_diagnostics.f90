! =============================================================================
!                         Test MPI field diagnostics
! =============================================================================
program test_mpi_field_diagnostics
    use constants, only : zero, one, f12
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    use fields
    use field_diagnostics
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, vcell, vcelli, ngrid
    use timer
    implicit none

    logical :: passed = .true.

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    call register_timer('field stats', field_stats_timer)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call update_parameters

    ! calls mpi_layout_init internally
    call field_alloc

    volg = vcell + one
    nparg = 1
    nsparg = 1
    velog = f12

    call calculate_field_diagnostics

    if (mpi_rank == mpi_master) then
        passed = (passed .and. (dabs(field_stats(IDX_RMS_V) - vcelli) == zero))
        passed = (passed .and. (dabs(field_stats(IDX_ABSERR_V) - vcelli) == zero))
        passed = (passed .and. (field_stats(IDX_MAX_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_MIN_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_AVG_NPAR) == one))
        passed = (passed .and. (field_stats(IDX_AVG_NSPAR) == one))
        passed = (passed .and. (field_stats(IDX_KEG) == 0.375d0 * dble(ngrid) * (vcell + one)))
    endif

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI field diagnostics', passed)
    endif

end program test_mpi_field_diagnostics
