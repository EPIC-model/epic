! =============================================================================
!                       Test ellipsoid split
!
!         This unit test checks the splitting of an ellipsoid. The ellipsoid
!         is centred at the origin with an orientation of 45 degrees.
! =============================================================================
program test_ellipsoid_split
    use unit_test
    use constants, only : pi, zero, one, three, four, five, ten, f12, f14, f34
    use parcel_container
    use parcel_split_mod, only : parcel_split, split_timer
    use parameters, only : update_parameters, nx, ny, nz, extent, lower, vmax
    use timer
    implicit none

    double precision, parameter :: lam = five
    double precision, parameter :: theta = f14 * pi, phi = pi / three !f14 * pi, phi = f14 * pi
    double precision, parameter :: st = dsin(theta), ct = dcos(theta)
    double precision, parameter :: sp = dsin(phi), cp = dcos(phi)
    double precision :: B11, B12, B13, B22, B23, B33, abc, ab, evec(3)
    double precision :: h, pos(2, 3), error, a2, b2, c2

    nx = 10
    ny = 10
    nz = 10
    extent = (/ten, ten, ten/)
    lower = (/-five, -five, -five/)
    call update_parameters

    vmax = one

    call register_timer('parcel split', split_timer)

    n_parcels = 1
    call parcel_alloc(2)

    abc = one
    ab = f34 * one

    a2 = ab * lam
    b2 = ab / lam
    c2 = (abc / ab) ** 2

    print *, "exact:", a2, b2, c2

    print *, "abc:", abc

    parcels%position(1, :) = zero
    parcels%volume(1) = four / three * abc * pi
    parcels%buoyancy(1) = one
#ifndef ENABLE_DRY_MODE
    parcels%humidity(1) = one
#endif

    B11 = a2 * st ** 2 * cp ** 2 + b2 * ct ** 2 * cp ** 2 + c2 * sp ** 2
    B12 = a2 * st ** 2 * sp * cp + b2 * ct ** 2 * sp * cp + c2 * sp * cp
    B13 = a2 * st * ct * cp - b2 * st * ct * cp
    B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 * sp ** 2 + c2 * cp ** 2
    B23 = a2 * st * ct * sp - b2 * st * ct * sp
    B33 = a2 * ct ** 2 + b2 * st ** 2

    print *, "B11:", B11
    print *, "B12:", B12
    print *, "B13:", B13
    print *, "B22:", B22
    print *, "B23:", B23
    print *, "B33:", B33


    parcels%B(1, 1) = B11
    parcels%B(1, 2) = B12
    parcels%B(1, 3) = B13
    parcels%B(1, 4) = B22
    parcels%B(1, 5) = B23

    ! analytic split
    h = f12 * dsqrt(three / five * a2)
    B11 = B11 - f34 * a2 * st ** 2 * cp ** 2
    B12 = B12 - f34 * a2 * st ** 2 * sp * cp
    B13 = B13 - f34 * a2 * st * ct * cp
    B22 = B22 - f34 * a2 * st ** 2 * sp ** 2
    B23 = B23 - f34 * a2 * st * ct * sp

    error = zero

    if ((a2 > b2) .and. (a2 > c2)) then
        evec = (/st * cp, st * sp, ct/)
    else if ((b2 > a2) .and. (b2 > c2)) then
        evec = (/ct * cp, ct * sp, -st/)
    else if ((c2 > a2) .and. (c2 > b2)) then
        evec = (/-sp, cp, zero/)
    else
        error = one
    endif

    pos(1, :) = parcels%position(1, :) + h * evec
    pos(2, :) = parcels%position(1, :) - h * evec

    ! numerical split
    call parcel_split(parcels, threshold=four)

    !
    ! check result
    !


    ! first parcel
    error = max(error, abs(parcels%B(1, 1) - B11))
    print *, "error", parcels%B(1, 1), B11
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, abs(parcels%B(1, 3) - B13))
    error = max(error, abs(parcels%B(1, 4) - B22))
    error = max(error, abs(parcels%B(1, 5) - B23))
    error = max(error, sum(abs(pos(1, :) - parcels%position(1, :))))
    print *, "error", error
    error = max(error, abs(f12 * four / three * abc * pi - parcels%volume(1)))
    error = max(error, abs(parcels%buoyancy(1) - one))
#ifndef ENABLE_DRY_MODE
    error = max(error, abs(parcels%humidity(1) - one))
#endif

    ! second parcel
    error = max(error, abs(parcels%B(2, 1) - B11))
    error = max(error, abs(parcels%B(2, 2) - B12))
    error = max(error, abs(parcels%B(2, 3) - B13))
    error = max(error, abs(parcels%B(2, 4) - B22))
    error = max(error, abs(parcels%B(2, 5) - B23))
    error = max(error, sum(abs(pos(2, :) - parcels%position(2, :))))
    error = max(error, abs(f12 * four / three * abc * pi - parcels%volume(2)))
    error = max(error, dble(abs(n_parcels - 2)))
    error = max(error, abs(parcels%buoyancy(2) - one))
#ifndef ENABLE_DRY_MODE
    error = max(error, abs(parcels%humidity(2) - one))
#endif
    call print_result_dp('Test ellipsoid split', error)

    call parcel_dealloc

end program test_ellipsoid_split
