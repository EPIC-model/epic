! =============================================================================
!                         Test MPI parcel diagnostics
! =============================================================================
program test_mpi_parcel_diagnostics
    use constants, only : zero, one, f12, f23
    use unit_test
    use mpi_environment
    use mpi_layout
    use parcels_mod, only : parcels
    use parcel_init, only : parcel_default, init_timer
    use parcel_diagnostics
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, vcell, set_vmin
    use mpi_timer
    implicit none

    logical                       :: passed = .true.
    integer, parameter            :: n_per_dim = 2
    double precision              :: total_vol

    call mpi_env_initialise

    passed = (world%err == 0)

    call register_timer('parcel stats', parcel_stats_timer)
    call register_timer('parcel init', init_timer)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    ! set to make all parcels smaller than vmin
    call set_vmin(vcell)

    call parcel_default

    parcels%volume(1:parcels%local_num) = vcell / dble(n_per_dim ** 3)
    parcels%B(:, 1:parcels%local_num) = zero
    parcels%B(1, 1:parcels%local_num) = parcels%get_abc(parcels%volume(1:parcels%local_num)) ** f23
    parcels%B(4, 1:parcels%local_num) = parcels%B(1, 1:parcels%local_num)
    parcels%vorticity(:, 1:parcels%local_num) = f12
    parcels%delta_pos(:, 1:parcels%local_num)  = f12
    parcels%buoyancy(1:parcels%local_num) = parcels%position(3, 1:parcels%local_num)

    call calculate_parcel_diagnostics

    if (world%rank == world%root) then
        total_vol = dble(parcels%total_num) * parcels%volume(1)
        passed = (passed .and. (dabs(parcel_stats(IDX_KE) - 0.375d0 * total_vol) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_N_SMALL) - parcels%total_num) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_AVG_LAM) - one) < 1.0e-13))
        passed = (passed .and. (parcel_stats(IDX_STD_LAM) < 1.0e-7))
        passed = (passed .and. (parcel_stats(IDX_STD_VOL) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_SUM_VOL) - total_vol) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_NTOT_PAR) - parcels%total_num) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_XI) - f12) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ETA) - f12) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ZETA) - f12) < 1.0e-15))
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel diagnostics', passed)
    endif

end program test_mpi_parcel_diagnostics
