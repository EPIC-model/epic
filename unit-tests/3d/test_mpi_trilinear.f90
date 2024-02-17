! =============================================================================
!                       Test trilinear interpolation
!
!         This unit test checks the trilinear interpolation par2grid
! =============================================================================
program test_mpi_trilinear
    use unit_test
    use options, only : parcel
    use mpi_environment
    use constants, only : pi, zero, one, f12, f23, f32
    use parcel_init, only : parcel_default, init_timer
    use mpi_layout
    use parcels_mod, only : parcels, top_parcels, bot_parcels
    use parcel_interpl, only : par2grid, par2grid_timer, halo_swap_timer
    use parameters, only : lower, update_parameters, nx, ny, nz, ngrid, extent
    use fields, only : vortg, field_alloc
    use field_ops, only : get_sum
    use mpi_timer
    implicit none

    double precision :: error
    logical          :: passed = .true.

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-f32, -f32, -f32/)
    extent =  (/0.4d0, 0.4d0, 0.4d0/)

    call register_timer('par2grid', par2grid_timer)
    call register_timer('halo swap', halo_swap_timer)
    call register_timer('parcel init', init_timer)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    parcel%min_vratio = 27.0d0
    parcel%n_per_cell = 27
    parcel%n_surf_per_cell = 9

    call update_parameters

    call field_alloc

    call parcel_default

    parcels%vorticity(1, 1:parcels%local_num) = 1.5d0
    bot_parcels%vorticity(1, 1:bot_parcels%local_num) = 1.5d0
    top_parcels%vorticity(1, 1:top_parcels%local_num) = 1.5d0

    call par2grid

    error = abs(get_sum(vortg(:, :, :, 1)) - 1.5d0 * dble(ngrid))

    passed = (passed .and. (error < dble(1.0e-15)))

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI trilinear (par2grid)', passed)
    endif

end program test_mpi_trilinear
