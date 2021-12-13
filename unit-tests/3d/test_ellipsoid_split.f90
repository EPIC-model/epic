! =============================================================================
!                       Test ellipsoid split
!
!         This unit test checks the splitting of an ellipsoid. The ellipsoid
!         is centred at the origin and the orientation is varied.
! =============================================================================
program test_ellipsoid_split
    use unit_test
    use constants, only : pi, zero, one, three, four, five, ten, f12, f14, f34
    use parcel_container
    use parcel_ellipsoid, only : get_B33
    use parcel_split_mod, only : parcel_split, split_timer
    use parameters, only : update_parameters, nx, ny, nz, extent, lower, vmax
    use timer
    implicit none

    double precision, parameter :: lam = five
    double precision :: theta, phi, st, ct, sp, cp
    double precision :: B11, B12, B13, B22, B23, B33, abc, ab, evec(3)
    double precision :: h, pos(3, 2), error, a2, b2, c2
    integer :: i, j

    nx = 10
    ny = 10
    nz = 10
    extent = (/ten, ten, ten/)
    lower = (/-five, -five, -five/)
    call update_parameters

    vmax = one

    call register_timer('parcel split', split_timer)

    call parcel_alloc(2)

    abc = one
    ab = f34

    a2 = ab * lam
    b2 = ab / lam
    c2 = (abc / ab) ** 2

    error = zero

    do i = 0, 7
        theta = dble(i) * f14 * pi
        st = dsin(theta)
        ct = dcos(theta)
        do j = 0, 7
            phi = dble(j) * f14 * pi
            sp = dsin(phi)
            cp = dcos(phi)

            call setup_parcels

            ! numerical split
            call parcel_split(parcels, threshold=four)

            call check_result
        enddo
    enddo

    call print_result_dp('Test ellipsoid split', error, atol=1.0e-14)

    call parcel_dealloc


    contains

        subroutine setup_parcels
            n_parcels = 1

            parcels%position(:, 1) = zero
            parcels%volume(1) = four / three * abc * pi
            parcels%buoyancy(1) = one
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = one
#endif
            ! 7 Nov 2021
            ! https://mathworld.wolfram.com/SphericalCoordinates.html
            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp
            B33 = a2 * cp ** 2 + c2 * sp ** 2

            parcels%B(1, 1) = B11
            parcels%B(2, 1) = B12
            parcels%B(3, 1) = B13
            parcels%B(4, 1) = B22
            parcels%B(5, 1) = B23
        end subroutine setup_parcels

        subroutine check_result
            ! analytic split
            h = f12 * dsqrt(three / five * a2)
            B11 = B11 - f34 * a2 * ct ** 2 * sp ** 2
            B12 = B12 - f34 * a2 * st * ct * sp ** 2
            B13 = B13 - f34 * a2 * ct * sp * cp
            B22 = B22 - f34 * a2 * st ** 2 * sp ** 2
            B23 = B23 - f34 * a2 * st * sp * cp
            B33 = B33 - f34 * a2 * cp ** 2

            if ((a2 > b2) .and. (a2 > c2)) then
                evec = (/ct * sp, st * sp, cp/)
            else if ((b2 > a2) .and. (b2 > c2)) then
                evec = (/-st, ct, zero/)
            else if ((c2 > a2) .and. (c2 > b2)) then
                evec = (/ct * cp, st * cp, -sp/)
            else
                error = one
            endif

            pos(:, 1) = h * evec
            pos(:, 2) = - h * evec

            ! exchange position
            if (sum(abs(pos(:, 1) - parcels%position(:, 1))) > 1.0e-13) then
                pos = -pos
            endif

            ! first parcel
            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(2, 1) - B12))
            error = max(error, abs(parcels%B(3, 1) - B13))
            error = max(error, abs(parcels%B(4, 1) - B22))
            error = max(error, abs(parcels%B(5, 1) - B23))
            error = max(error, abs(get_B33(parcels%B(:, 1), parcels%volume(1)) - B33))
            error = max(error, sum(abs(pos(:, 1) - parcels%position(:, 1))))
            error = max(error, abs(f12 * four / three * abc * pi - parcels%volume(1)))
            error = max(error, abs(parcels%buoyancy(1) - one))
#ifndef ENABLE_DRY_MODE
            error = max(error, abs(parcels%humidity(1) - one))
#endif

            ! second parcel
            error = max(error, abs(parcels%B(1, 2) - B11))
            error = max(error, abs(parcels%B(2, 2) - B12))
            error = max(error, abs(parcels%B(3, 2) - B13))
            error = max(error, abs(parcels%B(4, 2) - B22))
            error = max(error, abs(parcels%B(5, 2) - B23))
            error = max(error, abs(get_B33(parcels%B(:, 2), parcels%volume(2)) - B33))
            error = max(error, sum(abs(pos(:, 2) - parcels%position(:, 2))))
            error = max(error, abs(f12 * four / three * abc * pi - parcels%volume(2)))
            error = max(error, dble(abs(n_parcels - 2)))
            error = max(error, abs(parcels%buoyancy(2) - one))
#ifndef ENABLE_DRY_MODE
            error = max(error, abs(parcels%humidity(2) - one))
#endif
        end subroutine check_result

end program test_ellipsoid_split
