! =============================================================================
!                       Test 3D parcel correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_parcel_correction_3d
    use unit_test
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f23, f32, two, four, f12, f18
    use parcel_container
    use parcel_correction
    use parcel_bc
    use parcel_interpl, only : vol2grid
    use parcel_ellipsoid, only : get_abc
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, ny, nz, dx
    use fields, only : volg
    use timer
    implicit none

    double precision :: final_error, init_error, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    double precision :: q, m, delta

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)

    call parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)
    call register_timer('gradient correction', grad_corr_timer)
    call register_timer('vorticity correction', vort_corr_timer)


    nx = 32
    ny = 32
    nz = 32

    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

    n_parcels = 8*nx*ny*nz
    call parcel_alloc(n_parcels)

    parcel%n_per_cell = 8
    call init_regular_positions

    ! 0 --> -delta
    ! 1 -->  delta
    ! y = mx + q
    ! -delta = q
    ! delta = m - delta --> m = 2 * delta
    delta = f14 * dx(1)

    q = - delta
    m = two * delta

    ! add some deviation
    do n = 1, n_parcels
        call random_number(val)

        tmp = m * val + q
        parcels%position(1, n) = parcels%position(1, n) + tmp

        call random_number(val)
        tmp = m * val + q
        parcels%position(2, n) = parcels%position(2, n) + tmp

        call random_number(val)
        tmp = m * val + q
        parcels%position(3, n) = parcels%position(3, n) + tmp

        call apply_periodic_bc(parcels%position(:, n))
        call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
    enddo

    volg = zero

    parcels%volume = f18 * vcell

    parcels%B(:, :) = zero

    ! b11
    parcels%B(1, :) = get_abc(parcels%volume) ** f23

    ! b22
    parcels%B(4, :) = parcels%B(1, :)

    call vol2grid

    init_error = sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1))

    if (l_verbose) then
        open(unit=12, file='test_parcel_correction.asc', status='replace')
        write(12,*) '# iteration, average error, max absolute error'
        write(12,*) 0, init_error, maxval(abs(volg(0:nz, :, :) / vcell - one))
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        call vol2grid
        call apply_gradient(1.80d0,0.5d0)
        if (l_verbose) then
            call vol2grid
            write(12,*) i, sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1)), &
                           maxval(abs(volg(0:nz, :, :) / vcell - one))
        endif
    enddo

    if (l_verbose) then
        close(12)
    endif

    call vol2grid

    final_error = sum(abs(volg(0:nz, :, :) / vcell - one)) / (nx * ny * (nz+1))

    call print_result_dp('Test laplace and gradient 3D', final_error, init_error)

    deallocate(volg)

end program test_parcel_correction_3d

