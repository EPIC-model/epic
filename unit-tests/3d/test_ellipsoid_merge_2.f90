! =============================================================================
!                       Test ellipsoid multi merge
!
!         This unit test checks the merging of two ellipsoids. It is the
!         reverse test of test_ellipsoid_split.f90.
! =============================================================================
program test_ellipsoid_merge_2
    use unit_test
    use constants, only : pi, one, f12, f14, f13, f23, f34, two, four, five, ten
    use parcel_container
    use parcel_merge, only : merge_parcels, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, nx, ny, nz, lower, extent
    use parcel_ellipsoid
    use parcel_nearest
    use timer
    implicit none

    double precision :: a2, b2, c2, abc, ab, error, lam
    double precision :: theta, phi, st, ct, sp, cp
    integer          :: i, j

    nx = 1
    ny = 1
    nz = 1
    lower  = (/-pi / two, -pi /two, -pi /two/)
    extent = (/pi, pi, pi/)

    call register_timer('parcel merge', merge_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters

    call parcel_alloc(2)

    abc = one
    ab = f34

    lam = three
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

            call merge_parcels(parcels)

            call check_result
        enddo
    enddo

    call print_result_dp('Test ellipsoid merge 2', error, atol=dble(2.0e-15))

    call parcel_dealloc

    contains

        subroutine setup_parcels
            double precision :: h, B11, B12, B13, B22, B23, B33, evec(3)
            n_parcels = 2

            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp
            B33 = a2 * cp ** 2 + c2 * sp ** 2

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
            endif

            ! first parcel
            parcels%position(1, :) = h * evec
            parcels%volume(1) = f12 * dsqrt(a2 * b2 * c2) * four * pi / three
            parcels%B(1, 1) = B11
            parcels%B(1, 2) = B12
            parcels%B(1, 3) = B13
            parcels%B(1, 4) = B22
            parcels%B(1, 5) = B23
            parcels%buoyancy(1) = 1.8d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = 1.2d0
#endif
            ! second parcel
            parcels%position(2, :) = - h * evec
            parcels%volume(2) = parcels%volume(1)
            parcels%B(2, :) = parcels%B(1, :)
            parcels%buoyancy(2) = 1.4d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(2) = 1.1d0
#endif
        end subroutine setup_parcels

        subroutine check_result
            double precision :: B11, B12, B13, B22, B23, B33, vol
            double precision :: hum, buoy

            ! reference solution
            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp
            B33 = a2 * cp ** 2 + c2 * sp ** 2

            abc = dsqrt(a2 * b2 * c2)
            vol = abc * four * pi / three

            buoy = f12 * (1.8d0 + 1.4d0)
            hum  = f12 * (1.2d0 + 1.1d0)

            error = max(error, abs(dble(n_parcels - 1)))
            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(1, 2) - B12))
            error = max(error, abs(parcels%B(1, 3) - B13))
            error = max(error, abs(parcels%B(1, 4) - B22))
            error = max(error, abs(parcels%B(1, 5) - B23))
            error = max(error, abs(get_B33(parcels%B(1, :), &
                                           parcels%volume(1)) - B33))
            error = max(error, sum(abs(parcels%position(1, :))))
            error = max(error, abs(parcels%volume(1) - vol))
            error = max(error, abs(parcels%buoyancy(1) - buoy))
#ifndef ENABLE_DRY_MODE
            error = max(error, abs(parcels%humidity(1) - hum))
#endif
        end subroutine check_result
end program test_ellipsoid_merge_2
