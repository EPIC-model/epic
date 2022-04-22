! =============================================================================
!                         Test MPI field diagnostics
! =============================================================================
program test_field_diagnostics
    use constants, only : zero, one
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
    velog = one

    call calculate_field_diagnostics

    if (mpi_rank == mpi_master) then
        print *, ngrid, (box%hi(3)-box%lo(3) +1) * (box%hi(2)-box%lo(2)+1) * (box%hi(1)-box%lo(1)+1)
        print *, vcelli, field_stats(IDX_RMS_V)
!         passed = (passed .and. (field_stats(IDX_RMS_V) == one)
!         passed = (passed .and. (field_stats(IDX_ABSERR_V) == one)
!         passed = (passed .and. (field_stats(IDX_MAX_NPAR) == one)
!         passed = (passed .and. (field_stats(IDX_MIN_NPAR) == one)
!         passed = (passed .and. (field_stats(IDX_AVG_NPAR) == one)
!         passed = (passed .and. (field_stats(IDX_AVG_NSPAR) == one)
!         passed = (passed .and. (field_stats(IDX_KEG) == one)

    endif

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI field diagnostics', passed)
    endif

end program test_field_diagnostics
