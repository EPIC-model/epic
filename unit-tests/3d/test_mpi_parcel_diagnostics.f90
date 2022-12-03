! =============================================================================
!                         Test MPI parcel diagnostics
! =============================================================================
program test_mpi_parcel_diagnostics
    use constants, only : zero, one, f12, f23
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    use fields, only : field_alloc
    use parcel_container
    use parcel_diagnostics
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, vcell, dx, set_vmin
    use timer
    implicit none

    logical                       :: passed = .true.
    integer, parameter            :: n_per_dim = 2
    integer                       :: ix, iy, iz, i, j, k, l, n_total
    double precision              :: im, corner(3), total_vol
    double precision, allocatable :: velocity(:, :)

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    call register_timer('parcel stats', parcel_stats_timer)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call update_parameters

    ! set to make all parcels smaller than vmin
    call set_vmin(vcell)

    ! calls mpi_layout_init internally
    call field_alloc

    n_parcels = n_per_dim ** 3 * nz * (box%hi(2)-box%lo(2)+1) * (box%hi(1)-box%lo(1)+1)
    n_total = n_per_dim ** 3 * nz * ny * nx
    call parcel_alloc(n_parcels)

    allocate(velocity(3, n_parcels))

    im = one / dble(n_per_dim)

    l = 1
    do iz = 0, nz-1
        do iy = box%lo(2), box%hi(2)
            do ix = box%lo(1), box%hi(1)
                corner = lower + dble((/ix, iy, iz/)) * dx
                do k = 1, n_per_dim
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                            parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                            parcels%buoyancy(l) = parcels%position(3, l)
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    parcels%volume(1:n_parcels) = vcell / dble(n_per_dim ** 3)
    parcels%B(:, 1:n_parcels) = zero
    parcels%B(1, 1:n_parcels) = get_abc(parcels%volume(1:n_parcels)) ** f23
    parcels%B(4, 1:n_parcels) = parcels%B(1, 1:n_parcels)
    parcels%vorticity(:, 1:n_parcels) = f12
    velocity(:, 1:n_parcels)  = f12

    ! calculates reference potential energy
    call calculate_peref

    call calculate_parcel_diagnostics(velocity)

    if (mpi_rank == mpi_master) then
        total_vol = dble(n_total) * parcels%volume(1)
        passed = (passed .and. (dabs(parcel_stats(IDX_KE) - 0.375d0 * total_vol) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_N_SMALL) - n_total) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_AVG_LAM) - one) < 1.0e-13))
        passed = (passed .and. (parcel_stats(IDX_STD_LAM) < 1.0e-7))
        passed = (passed .and. (parcel_stats(IDX_STD_VOL) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_SUM_VOL) - total_vol) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_NTOT_PAR) - n_total) == zero))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_XI) - f12) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ETA) - f12) < 1.0e-15))
        passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ZETA) - f12) < 1.0e-15))
    endif

    deallocate(velocity)

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI parcel diagnostics', passed)
    endif

end program test_mpi_parcel_diagnostics
