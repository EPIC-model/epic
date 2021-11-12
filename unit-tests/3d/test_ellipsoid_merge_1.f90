! =============================================================================
!                       Test ellipsoid merge
!
!         This unit test checks the merging of 9 ellipsoids. The biggest
!         ellipsoid is located at the origin. The smaller ellipsoids are located
!         on all sides. The final ellipsoid is a circle located at the
!         origin.
! =============================================================================
program test_ellipsoid_merge_1
    use unit_test
    use constants, only : pi, three, four, five, eight, ten, f23
    use parcel_container
    use parcel_merge, only : merge_parcels, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz
    use parcel_ellipsoid
    use parcel_nearest
    use timer
    implicit none

    double precision :: a1b1c1, a2b2c2, error

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

    call parcel_alloc(9)

    a1b1c1 = 0.49d0
    a2b2c2 = 0.25d0

    call parcel_setup

    call merge_parcels(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipsoid merge 1', error, atol=dble(2.0e-15))

    call parcel_dealloc

    contains

        subroutine parcel_setup
            integer :: i

            n_parcels = 9
            parcels%position(1, :) = zero
            parcels%volume(1) = a1b1c1 * four * pi / three
            parcels%B(1, :) = zero
            parcels%B(1, 1) = a1b1c1 ** f23
            parcels%B(1, 4) = a1b1c1 ** f23
            parcels%buoyancy(1) = 1.5d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = 1.3d0
#endif

            do i = 0, 1
                ! small parcel below front/back left
                parcels%position(2+4*i, 1) = -0.6d0
                parcels%position(2+4*i, 2) = -0.6d0 + 1.2d0 * dble(i)
                parcels%position(2+4*i, 3) = -0.6d0
                parcels%volume(2+4*i) = a2b2c2 * four * pi / three
                parcels%B(2+4*i, :) = zero
                parcels%B(2+4*i, 1) = a2b2c2 ** f23
                parcels%B(2+4*i, 4) = a2b2c2 ** f23
                parcels%buoyancy(2+4*i) = 1.8d0
#ifndef ENABLE_DRY_MODE
                parcels%humidity(2+4*i) = 1.2d0
#endif
                ! small parcel below front/back right
                parcels%position(3+4*i, 1) =  0.6d0
                parcels%position(3+4*i, 2) = -0.6d0 + 1.2d0 * dble(i)
                parcels%position(3+4*i, 3) = -0.6d0
                parcels%volume(3+4*i) = a2b2c2 * four * pi / three
                parcels%B(3+4*i, :) = zero
                parcels%B(3+4*i, 1) = a2b2c2 ** f23
                parcels%B(3+4*i, 4) = a2b2c2 ** f23
                parcels%buoyancy(3+4*i) = 1.4d0
#ifndef ENABLE_DRY_MODE
                parcels%humidity(3+4*i) = 1.1d0
#endif
                ! small parcel above front/back left
                parcels%position(4+4*i, 1) = -0.6d0
                parcels%position(4+4*i, 2) = -0.6d0 + 1.2d0 * dble(i)
                parcels%position(4+4*i, 3) =  0.6d0
                parcels%volume(4+4*i) = a2b2c2 * four * pi / three
                parcels%B(4+4*i, :) = zero
                parcels%B(4+4*i, 1) = a2b2c2 ** f23
                parcels%B(4+4*i, 4) = a2b2c2 ** f23
                parcels%buoyancy(4+4*i) = 1.7d0
#ifndef ENABLE_DRY_MODE
                parcels%humidity(4+4*i) = 1.0d0
#endif
                ! small parcel above front/back right
                parcels%position(5+4*i, 1) =  0.6d0
                parcels%position(5+4*i, 2) = -0.6d0 + 1.2d0 * dble(i)
                parcels%position(5+4*i, 3) =  0.6d0
                parcels%volume(5+4*i) = a2b2c2 * four * pi / three
                parcels%B(5+4*i, :) = zero
                parcels%B(5+4*i, 1) = a2b2c2 ** f23
                parcels%B(5+4*i, 4) = a2b2c2 ** f23
                parcels%buoyancy(5+4*i) = 1.5d0
#ifndef ENABLE_DRY_MODE
                parcels%humidity(5+4*i) = 1.4d0
#endif
            enddo
        end subroutine parcel_setup

        function eval_max_error() result(max_err)
            double precision :: abc, B11, B12, B13, B22, B23, B33, vol, buoy
#ifndef ENABLE_DRY_MODE
            double precision :: hum
#endif
            double precision :: max_err

            ! reference solution
            abc = a1b1c1 + eight * a2b2c2  ! a == b == c since it is a sphere
            B11 = abc ** f23
            B12 = zero
            B13 = zero
            B22 = abc ** f23
            B23 = zero
            B33 = abc ** f23
            vol = abc * four * pi / three

            buoy = (1.5d0 * a1b1c1 + two * (1.8d0 + 1.4d0 + 1.7d0 + 1.5d0) * a2b2c2) / abc
#ifndef ENABLE_DRY_MODE
            hum  = (1.3d0 * a1b1c1 + two * (1.2d0 + 1.1d0 + 1.0d0 + 1.4d0) * a2b2c2) / abc
#endif
            max_err = zero

            max_err = max(max_err, abs(dble(n_parcels - 1)))
            max_err = max(max_err, abs(parcels%B(1, 1) - B11))
            max_err = max(max_err, abs(parcels%B(1, 2) - B12))
            max_err = max(max_err, abs(parcels%B(1, 3) - B13))
            max_err = max(max_err, abs(parcels%B(1, 4) - B22))
            max_err = max(max_err, abs(parcels%B(1, 5) - B23))
            max_err = max(max_err, abs(get_B33(parcels%B(1, :), &
                                               parcels%volume(1)) - B33))
            max_err = max(max_err, sum(abs(parcels%position(1, :))))
            max_err = max(max_err, abs(parcels%volume(1) - vol))
            max_err = max(max_err, abs(parcels%buoyancy(1) - buoy))
#ifndef ENABLE_DRY_MODE
            max_err = max(max_err, abs(parcels%humidity(1) - hum))
#endif
        end function eval_max_error



end program test_ellipsoid_merge_1
