! =============================================================================
!                       Test grid2par interpolation
!
!                       This unit test checks grid2par
! =============================================================================
program test_mpi_grid2par
    use unit_test
    use options, only : parcel
    use mpi_environment
    use mpi_collectives
    use mpi_layout
    use constants, only : pi, zero, one, two, three, four, five, f12, f23
    use parcels_mod, only : parcels, bot_parcels, top_parcels
    use parcel_init, only : parcel_default, init_timer
    use parcel_interpl, only : grid2par, grid2par_timer
    use parameters, only : lower, update_parameters, nx, ny, nz, extent
    use fields, only : velog, vortg, velgradg, field_alloc
    use mpi_timer
    implicit none

    double precision :: error
    logical          :: passed = .true.
    integer          :: l

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-one, -one, -one/)
    extent =  (/two, two, two/)

    call register_timer('grid2par', grid2par_timer)
    call register_timer('parcel init', init_timer)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    parcel%min_vratio = 27.0d0
    parcel%n_per_cell = 27
    parcel%n_surf_per_cell = 9

    call update_parameters

    call field_alloc

    call parcel_default

    velog(:, :, :, 1) = one
    velog(:, :, :, 2) = two
    velog(:, :, :, 3) = three

    vortg(:, :, :, 1) = one
    vortg(:, :, :, 2) = two
    vortg(:, :, :, 3) = three

    do l = 1, 5
        velgradg(:, :, :, l) = dble(l)
    enddo

    ! we cannot check parcels%delta_vor since parcels%vorticity is used in grid2par
    call grid2par

    error = zero

    ! interior parcels:
    do l = 1, 3
        error = max(error, maxval(abs(parcels%delta_pos(l, 1:parcels%local_num) - dble(l))))
    enddo

    do l = 1, 5
        error = max(error, maxval(abs(parcels%strain(l, 1:parcels%local_num) - dble(l))))
    enddo


    ! surface parcels:
    do l = 1, 2
        error = max(error, maxval(dabs(bot_parcels%delta_pos(l, 1:bot_parcels%local_num) - dble(l))))
        error = max(error, maxval(dabs(top_parcels%delta_pos(l, 1:top_parcels%local_num) - dble(l))))
    enddo

    ! du/dx = 1 and du/dy = 2
    do l = 1, 2
        error = max(error, maxval(dabs(bot_parcels%strain(l, 1:bot_parcels%local_num) - dble(l))))
        error = max(error, maxval(dabs(top_parcels%strain(l, 1:top_parcels%local_num) - dble(l))))
    enddo

    ! dv/dx = \zeta + du/dy = 3 + 2 = 5
    error = max(error, maxval(dabs(bot_parcels%strain(3, 1:bot_parcels%local_num) - 5.0d0)))
    error = max(error, maxval(dabs(top_parcels%strain(3, 1:top_parcels%local_num) - 5.0d0)))

    ! dv/dy = 3
    error = max(error, maxval(dabs(bot_parcels%strain(4, 1:bot_parcels%local_num) - 3.0d0)))
    error = max(error, maxval(dabs(top_parcels%strain(4, 1:top_parcels%local_num) - 3.0d0)))

    call mpi_blocking_reduce(error, MPI_MAX, world)

    passed = (passed .and. (error < dble(1.0e-15)))

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI grid2par', passed)
    endif

end program test_mpi_grid2par
