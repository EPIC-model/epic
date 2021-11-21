! =============================================================================
!                       Test grid2par interpolation
!
!                       This unit test checks grid2par
! =============================================================================
program test_trilinear
    use unit_test
    use constants, only : pi, zero, one, two, three, f12, f23
    use parcel_container
    use parcel_interpl, only : grid2par, grid2par_timer
    use parcel_ellipsoid, only : get_abc
    use parameters, only : lower, update_parameters, vcell, dx, nx, ny, nz, ngrid
    use fields, only : velog, velgradg, vtend, field_alloc
    use timer
    implicit none

    double precision              :: error
    integer                       :: ix, iy, iz, i, j, k, l, n_per_dim
    double precision              :: im, corner(3)
    double precision, allocatable :: vel(:, :), vor(:, :), vgrad(:, :)

    nx = 32
    ny = 32
    nz = 32
    lower  = (/-one, -one, -one/)
    extent =  (/two, two, two/)

    call register_timer('grid2par', grid2par_timer)

    call update_parameters

    call field_alloc

    n_per_dim = 3

    n_parcels = n_per_dim ** 3 * nx *ny *nz
    call parcel_alloc(n_parcels)

    allocate(vel(n_parcels, 3))
    allocate(vor(n_parcels, 3))
    allocate(vgrad(n_parcels, 9))


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

    vtend(:, :, :, 1) = one
    vtend(:, :, :, 2) = two
    vtend(:, :, :, 3) = three

    velog(:, :, :, 1) = one
    velog(:, :, :, 2) = two
    velog(:, :, :, 3) = three

    do l = 1, 5
        velgradg(:, :, :, l) = dble(l)
    enddo

    call grid2par(vel, vor, vgrad)

    error = zero

    do l = 1, 3
        error = max(error, maxval(dabs(vel(1:n_parcels, l) - dble(l))))
        error = max(error, maxval(dabs(vor(1:n_parcels, l) - dble(l))))
    enddo

    do l = 1, 9
        error = max(error, maxval(dabs(vgrad(1:n_parcels, l) - dble(l))))
    enddo

    call print_result_dp('Test grid2par', error, atol=dble(1.0e-14))

    deallocate(vel)
    deallocate(vor)
    deallocate(vgrad)

end program test_trilinear
