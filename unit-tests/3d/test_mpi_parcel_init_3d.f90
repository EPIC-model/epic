! =============================================================================
!                       Test parcel initialisation in 3D
!
!         This unit test checks the parcel initialisation from fields.
! =============================================================================
program test_mpi_parcel_init_3d
    use unit_test
    use mpi_environment
    use options, only : parcel
    use field_mpi
    use constants, only : pi, zero, one, two, four, five, f12, f13, f23, f32
    use parcel_init, only : init_timer, parcel_default, init_parcels_from_grids
    use parcel_interpl, only : par2grid, par2grid_timer, halo_swap_timer
    use fields, only : tbuoyg, field_default
    use field_ops, only : get_rms, get_abs_max
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_layout, only : box, mpi_layout_init
    use mpi_timer
    implicit none

    double precision :: xg, yg, zg, facx, facy, facz, argx, argy, argz
    integer          :: ix, iy, iz
    double precision :: rms, rmserr, error
    logical          :: passed = .true.
    double precision, allocatable :: workg(:, :, :)
    double precision :: tol = 2.0d-2

    call mpi_env_initialise

    call parse_command_line

    passed = (world%err == 0)

    nx = 128
    ny = 128
    nz = 64
    lower = (/-four, -four, -two/)
    extent = (/8.0d0, 8.0d0, four/)

    call register_timer('parcel init', init_timer)
    call register_timer('par2grid', par2grid_timer)
    call register_timer('halo swap', halo_swap_timer)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    parcel%n_per_cell = 8

    call update_parameters

    call field_default

    !--------------------------------------------------------
    ! Define a gridded field "tbuoyg" (this can be arbitrary):
    facx = two * pi / extent(1)
    facy = two * pi / extent(2)
    facz = one / extent(3)
    do ix = box%lo(1), box%hi(1)
        xg = dx(1) * dble(ix)
        argx = facx * xg
        do iy = box%lo(2), box%hi(2)
            yg = dx(2) * dble(iy)
            argy = facy * yg
            do iz = 0, nz
                zg = dx(3) * dble(iz)
                argz = facz * zg
                tbuoyg(iz, iy, ix) = exp(argz) * sin(argx - argy + one) / &
                                     ((one + f12 * cos(argx)) * (one + f13 * cos(argy)))
            enddo
        enddo
    enddo

    call field_halo_swap_scalar(tbuoyg)

    !---------------------------
    ! Generate parcel attribute:
    call parcel_default

    call init_parcels_from_grids

    !
    ! check result
    !
    error = zero

    ! Copy buoyg since it is overwritten in par2grid after:
    allocate(workg(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    workg = tbuoyg

    call par2grid

    ! Compute max and rms errors:

    ! Rms value of original field
    rms = get_rms(workg)

    workg = tbuoyg - workg

    ! Rms error in reconstruction
    rmserr = get_rms(workg)

    error = max(error, rmserr)

    ! Relative rms error
    error = max(error, rmserr / rms)

    ! Maximum error
    error = max(error, get_abs_max(workg))

    passed = (passed .and. (error < tol))

    deallocate(workg)

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test parcel initialisation 3D', passed)
    endif

end program test_mpi_parcel_init_3d
