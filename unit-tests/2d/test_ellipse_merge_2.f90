! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the merging of three ellipses. The biggest
!         ellipse is located at the origin. The smaller ellipses are located
!         tangentially at 45 and 225 degrees. The final ellipse is an ellipse
!         located at the origin.
! =============================================================================
program test_ellipse_multi_merge_2
    use unit_test
    use constants, only : pi, one, two, four, five, ten
    use parcel_container
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, nx, nz, lower, extent
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

    call parcel_alloc(3)

    a1b1 = 1.44d0
    a2b2 = 0.25d0

    call parcel_setup

    call merge_ellipses(parcels)

    ! check result
    call eval_max_error

    call print_result_dp('Test ellipse group-merge 2', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            double precision :: d

            d = (sqrt(a1b1) + sqrt(a2b2)) * f12 * sqrt(two)

            n_parcels = 3
            parcels%position(:, 1) = zero
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(2, 1) = zero
            parcels%vorticity(1) = zero
            parcels%buoyancy(1) = 1.5d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = 1.3d0
#endif
            ! small parcel left
            parcels%position(1, 2) = -d
            parcels%position(2, 2) = -d
            parcels%volume(2) = a2b2 * pi
            parcels%B(1, 2) = a2b2
            parcels%B(2, 2) = zero
            parcels%vorticity(2) = zero
            parcels%buoyancy(2) = 1.8d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(2) = 1.2d0
#endif
            ! small parcel right
            parcels%position(1, 3) = d
            parcels%position(2, 3) = d
            parcels%volume(3) = a2b2 * pi
            parcels%B(1, 3) = a2b2
            parcels%B(2, 3) = zero
            parcels%vorticity(3) = zero
            parcels%buoyancy(3) = 1.4d0
#ifndef ENABLE_DRY_MODE
            parcels%humidity(3) = 1.1d0
#endif
        end subroutine parcel_setup

        subroutine eval_max_error
            double precision :: ab, B11, B12, B22, vol, d, factor
            double precision :: hum, buoy

            ! reference solution
            d = (sqrt(a1b1) + sqrt(a2b2)) * f12 * sqrt(two)
            ab = a1b1 + two * a2b2
            vol = ab * pi

            buoy = (1.5d0 * a1b1 + (1.8d0 + 1.4d0) * a2b2) / ab
            hum  = (1.3d0 * a1b1 + (1.2d0 + 1.1d0) * a2b2) / ab

            B11 = a1b1 ** 2 / ab + two * a2b2 / ab * (four * d ** 2 + a2b2)
            B12 = two * a2b2 / ab * (four * d ** 2)
            B22 = B11

            factor = ab / sqrt(B11 * B22 - B12 ** 2)
            B11 = B11 * factor
            B12 = B12 * factor
            B22 = B22 * factor

            error = zero
            error = max(error, abs(dble(n_parcels - 1)))
            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(2, 1) - B12))
            error = max(error, abs(get_B22(parcels%B(1, 1), &
                                           parcels%B(2, 1), &
                                           parcels%volume(1)) - B22))
            error = max(error, sum(abs(parcels%position(:, 1))))
            error = max(error, abs(parcels%volume(1) - vol))
            error = max(error, abs(parcels%buoyancy(1) - buoy))
#ifndef ENABLE_DRY_MODE
            error = max(error, abs(parcels%humidity(1) - hum))
#endif
        end subroutine eval_max_error
end program test_ellipse_multi_merge_2
