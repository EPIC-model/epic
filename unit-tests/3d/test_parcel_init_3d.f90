! =============================================================================
!                       Test parcel initialisation in 3D
!
!         This unit test checks the parcel initialisation from fields.
! =============================================================================
program test_parcel_init_3d
    use unit_test
    use constants, only : pi, zero, one, two, four, f12, f13, f23, f32
    use parcel_container
    use parcel_init, only : gen_parcel_scalar_attr, unit_test_parcel_init_alloc, init_timer
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipsoid, only : get_abc
    use fields, only : tbuoyg, field_default
    use parameters, only : update_parameters, dx, ncell, nx, ny, nz, lower, vcell
    use timer
    implicit none

    double precision  :: xg, yg, zg, facx, facy, facz, argx, argy, argz, v0
    integer :: i, ix, iy, iz, mx, my, mz
    double precision :: rms, rmserr, error
    double precision, allocatable :: workg(:, :, :)
    double precision :: tol = 1.0d-9

     !Number of parcels per grid box = nbgx*nbgz
    integer, parameter :: nbgx = 2, nbgy = nbgx, nbgz = nbgx

    !Fractions of grid cell used for placing parcels:
    double precision, parameter :: dxf = one / dble(nbgx), &
                                   dyf = one / dble(nbgy), &
                                   dzf = one / dble(nbgz)

    nx = 64
    ny = 64
    nz = 32
    lower = (/-four, -four, -two/)
    extent = (/8.0d0, 8.0d0, four/)

    call register_timer('parcel init', init_timer)
    call register_timer('par2grid', par2grid_timer)

    call update_parameters

    call field_default

    !Maximum number of parcels:

    !Total number of parcels:
    n_parcels = nbgx * nbgy * nbgz * ncell
    call parcel_alloc(n_parcels)


    !--------------------------------------------------------
    ! Define a gridded field "tbuoyg" (this can be arbitrary):
    facx = two * pi / extent(1)
    facy = two * pi / extent(2)
    facz = one / extent(3)
    do ix = 0, nx-1
        xg = dx(1) * dble(ix)
        argx = facx * xg
        do iy = 0, ny-1
            yg = dx(2) * dble(iy)
            argy = facy * yg
            do iz = 0, nz
                zg = dx(3) * dble(iz)
                argz = facz * zg
                tbuoyg(iz, iy, ix) = dexp(argz) * dsin(argx - argy + one) / &
                                     ((one + f12 * dcos(argx)) * (one + f13 * dcos(argy)))
            enddo
        enddo
    enddo

    !---------------------------------------------------------
    !Initialise parcel volume positions and volume fractions:
    v0 = dxf * dyf * dzf * vcell
    i = 0
    do iz = 0, nz-1
        do iy = 0, ny-1
            do ix = 0, nx-1
                do mz = 1, nbgz
                    do my = 1, nbgy
                        do mx = 1, nbgx
                            i = i + 1
                            parcels%position(1, i) = lower(1) + dx(1) * (dble(ix) + dxf * (dble(mx) - f12))
                            parcels%position(2, i) = lower(2) + dx(2) * (dble(iy) + dyf * (dble(my) - f12))
                            parcels%position(3, i) = lower(3) + dx(3) * (dble(iz) + dzf * (dble(mz) - f12))
                            parcels%volume(i) = v0
                            parcels%B(:, i) = zero
                            parcels%B(1, i) = get_abc(v0) ** f23
                            parcels%B(4, i) = parcels%B(1, i)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo


    ! Prepare for "gen_parcel_scalar_attr"
    call unit_test_parcel_init_alloc

    !---------------------------
    ! Generate parcel attribute:
    call gen_parcel_scalar_attr(tbuoyg, tol, parcels%buoyancy)

    !
    ! check result
    !
    error = zero

    ! Copy buoyg since it is overwritten in par2grid after:
    allocate(workg(-1:nz+1, 0:ny-1, 0:nx-1))
    workg = tbuoyg

    call par2grid

    ! Compute max and rms errors:

    ! Rms value of original field
    rms = dsqrt((f12 * sum(workg(0, :, :) ** 2 + workg(nz, :, :) ** 2) + &
                       sum(workg(1:nz-1, :, :) ** 2)) / dble(ncell))


    workg = tbuoyg - workg

    ! Rms error in reconstruction
    rmserr = dsqrt((f12 * sum(workg(0, :, :) ** 2 + workg(nz, :, :) ** 2) + &
                          sum(workg(1:nz-1, :, :) ** 2)) / dble(ncell))

    error = max(error, rmserr)

    ! Relative rms error
    error = max(error, rmserr / rms)

    ! Maximum error
    error = max(error, maxval(dabs(workg(0:nz, :, :))))

    call print_result_dp('Test parcel initialisation 3D', error, atol=two * tol)

    call parcel_dealloc

    deallocate(workg)

end program test_parcel_init_3d
