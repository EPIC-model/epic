! =============================================================================
!                       Test par2grid
!
!         This unit test checks the par2grid
! =============================================================================
program test_par2grid
    use unit_test
    use options, only : parcel
    use constants, only : pi, zero, one, two, f12, f23, f18
    use parcel_container
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipsoid, only : get_abc
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz
    use fields, only : vortg, field_default
    use timer
    implicit none

    double precision              :: error, x, y, z
    integer                       :: ix, iy, iz, n
    double precision              :: k, l, m, fk2l2, alpha, lam, l23
    double precision              :: cosmz,  sinmz, sinkxly, coskxly
    double precision, allocatable :: vortg_ref(:, :, :, :)

    nx = 32
    ny = 32
    nz = 32
    lower  = f12*pi*(/-one, -one, -one/)
    extent =  (/pi, pi, pi/)

    call register_timer('par2grid', par2grid_timer)

    allocate(vortg_ref(-1:nz+1, 0:ny-1, 0:nx-1, 3))

    call update_parameters

    call field_default

    n_parcels = 8*nx*ny*nz
    call parcel_alloc(n_parcels)

    parcel%n_per_cell = 8
    call init_regular_positions

    k = two
    l = two
    m = one

    alpha = dsqrt(k ** 2 + l ** 2 + m ** 2)
    fk2l2 = one / dble(k ** 2 + l ** 2)

    ! aspect ratio: lam = a / c
    lam = maxval((/dx(1) / dx(2), dx(2) / dx(1),   &
                   dx(1) / dx(3), dx(3) / dx(1),   &
                   dx(2) / dx(3), dx(3) / dx(2)/))

    do n = 1, n_parcels
        x = parcels%position(1, n)
        y = parcels%position(2, n)
        z = parcels%position(3, n)

        cosmz = dcos(m * z)
        sinmz = dsin(m * z)
        sinkxly = dsin(k * x + l * y)
        coskxly = dcos(k * x + l * y)

        parcels%vorticity(1, n) = alpha * fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
        parcels%vorticity(2, n) = alpha * fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
        parcels%vorticity(3, n) = alpha * cosmz * coskxly

        parcels%volume(n) = f18 * vcell
        parcels%B(:, n) = zero

        l23 = (lam * get_abc(parcels%volume(n))) ** f23

        ! b11
        parcels%B(1, n) = l23

        ! b22
        parcels%B(4, n) = l23
    enddo

    do ix = 0, nx-1
        x = lower(1) + ix * dx(1)
        do iy = 0, ny-1
            y = lower(2) + iy * dx(2)
            do iz = -1, nz+1
                z = lower(3) + iz * dx(3)

                cosmz = dcos(m * z)
                sinmz = dsin(m * z)
                sinkxly = dsin(k * x + l * y)
                coskxly = dcos(k * x + l * y)

                ! exact vorticity
                vortg_ref(iz, iy, ix, 1) = alpha * fk2l2 * (k * m * sinmz - l * alpha * cosmz) * sinkxly
                vortg_ref(iz, iy, ix, 2) = alpha * fk2l2 * (l * m * sinmz + k * alpha * cosmz) * sinkxly
                vortg_ref(iz, iy, ix, 3) = alpha * cosmz * coskxly

            enddo
        enddo
    enddo

    call par2grid(.false.)

    print *, nz, minval(vortg_ref(nz, :, :, 3)), maxval(vortg_ref(nz, :, :, 3))
    print *, nz, minval(vortg(nz, :, :, 3)), maxval(vortg(nz, :, :, 3))

    vortg(0:nz, :, :, 1) = dabs(vortg(0:nz, :, :, 1) - vortg_ref(0:nz, :, :, 1))
    vortg(0:nz, :, :, 2) = dabs(vortg(0:nz, :, :, 2) - vortg_ref(0:nz, :, :, 2))
    vortg(0:nz, :, :, 3) = dabs(vortg(0:nz, :, :, 3) - vortg_ref(0:nz, :, :, 3))
    print *, "abs max error in xi:", maxval(vortg(0:nz, :, :, 1))
    print *, "abs max error in eta:", maxval(vortg(0:nz, :, :, 2))
    print *, "abs max error in zeta:", maxval(vortg(0:nz, :, :, 3))
    print *, "location for xi:", maxloc(vortg(0:nz, :, :, 1))
    print *, "location for eta:", maxloc(vortg(0:nz, :, :, 2))
    print *, "location for zeta:", maxloc(vortg(0:nz, :, :, 3))

    error = zero
!     call print_result_dp('Test par2grid', error, atol=dble(3.0e-14))

    deallocate(vortg_ref)

end program test_par2grid
