! =============================================================================
!                       Test trilinear interpolation
!
!         This unit test checks the trilinear interpolation par2grid
! =============================================================================
program test_trilinear
    use unit_test
    use constants, only : pi, zero, one, f12, f23, f32
    use parcel_container
    use parcel_interpl, only : par2grid, par2grid_timer
    use parcel_ellipsoid, only : get_abc
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : volg, field_alloc
    use timer
    implicit none

    double precision :: error
    integer          :: ix, iy, iz, i, j, k, l, n_per_dim
    double precision :: im, corner(3)

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-f32, -f32, -f32/)
    extent =  (/0.4d0, 0.4d0, 0.4d0/)

    call register_timer('par2grid', par2grid_timer)

    call update_parameters

    call field_alloc

    n_per_dim = 3

    n_parcels = n_per_dim ** 3 * nx *ny *nz
    call parcel_alloc(n_parcels)


    im = one / dble(n_per_dim)

    l = 1
    do iz = 0, nz-1
        do iy = 0, ny-1
            do ix = 0, nx-1
                corner = lower + dble((/ix, iy, iz/)) * dx
                do k = 1, n_per_dim
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            parcels%position(l, 1) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(l, 2) = corner(2) + dx(2) * (dble(j) - f12) * im
                            parcels%position(l, 3) = corner(3) + dx(3) * (dble(k) - f12) * im
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    parcels%volume = vcell / dble(n_per_dim ** 3)

    parcels%B(:, :) = zero

    ! b11
    parcels%B(:, 1) = get_abc(parcels%volume(1:n_parcels)) ** f23

    ! b22
    parcels%B(:, 4) = parcels%B(:, 1)

    call par2grid

    error = abs(sum(volg(0:nz, :, :)) - dble(ngrid) * vcell)

    call print_result_dp('Test trilinear (par2grid)', error, atol=dble(3.0e-14))

end program test_trilinear
