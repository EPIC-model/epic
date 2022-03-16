! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the merging of five ellipses. The biggest
!         ellipse is located at the origin. The smaller ellipses are located
!         on all four sides. The final ellipse is a circle located at the
!         origin.
! =============================================================================
program test_ellipse_multi_merge_1
    use unit_test
    use constants, only : pi, four, five, ten
    use parcel_container
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_ellipse
    use parcel_nearest
    use timer
    implicit none

    double precision :: a1b1, a2b2, error

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('parcel merge', merge_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters

    call parcel_alloc(5)

    a1b1 = 0.49d0
    a2b2 = 0.25d0

    call parcel_setup

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse group-merge 1', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            n_parcels = 5
            parcels%position(:, 1) = zero
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(2, 1) = zero
            parcels%buoyancy(1) = 1.5d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = 1.3d0
#endif
            ! small parcel left
            parcels%position(1, 2) = -0.6d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = a2b2 * pi
            parcels%B(1, 2) = a2b2
            parcels%B(2, 2) = zero
            parcels%buoyancy(2) = 1.8d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(2) = 1.2d0
#endif
            ! small parcel right
            parcels%position(1, 3) = 0.6d0
            parcels%position(2, 3) = zero
            parcels%volume(3) = a2b2 * pi
            parcels%B(1, 3) = a2b2
            parcels%B(2, 3) = zero
            parcels%buoyancy(3) = 1.4d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(3) = 1.1d0
#endif
            ! small parcel below
            parcels%position(1, 4) = zero
            parcels%position(2, 4) = -0.6d0
            parcels%volume(4) = a2b2 * pi
            parcels%B(1, 4) = a2b2
            parcels%B(2, 4) = zero
            parcels%buoyancy(4) = 1.7d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(4) = 1.0d0
#endif
            ! small parcel above
            parcels%position(1, 5) = zero
            parcels%position(2, 5) = 0.6d0
            parcels%volume(5) = a2b2 * pi
            parcels%B(1, 5) = a2b2
            parcels%B(2, 5) = zero
            parcels%buoyancy(5) = 1.5d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(5) = 1.4d0
#endif
        end subroutine parcel_setup

        function eval_max_error() result(max_err)
            double precision :: ab, B11, B12, B22, vol, buoy
#ifndef ENABLE_DRY_MODE
            double precision :: hum
#endif
            double precision :: max_err

            ! reference solution
            ab = a1b1 + four * a2b2  ! a == b since it is a circle
            B11 = ab
            B12 = zero
            B22 = ab
            vol = ab * pi

            buoy = (1.5d0 * a1b1 + (1.8d0 + 1.4d0 + 1.7d0 + 1.5d0) * a2b2) / ab
#ifndef ENABLE_DRY_MODE
            hum  = (1.3d0 * a1b1 + (1.2d0 + 1.1d0 + 1.0d0 + 1.4d0) * a2b2) / ab
#endif
            max_err = zero

            max_err = max(max_err, abs(dble(n_parcels - 1)))
            max_err = max(max_err, abs(parcels%B(1, 1) - B11))
            max_err = max(max_err, abs(parcels%B(2, 1) - B12))
            max_err = max(max_err, abs(get_B22(parcels%B(1, 1), &
                                               parcels%B(2, 1), &
                                               parcels%volume(1)) - B22))
            max_err = max(max_err, sum(abs(parcels%position(:, 1))))
            max_err = max(max_err, abs(parcels%volume(1) - vol))
            max_err = max(max_err, abs(parcels%buoyancy(1) - buoy))
#ifndef ENABLE_DRY_MODE
            max_err = max(max_err, abs(parcels%humidity(1) - hum))
#endif
        end function eval_max_error



end program test_ellipse_multi_merge_1
