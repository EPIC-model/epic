! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the merging of five ellipses. The biggest
!         ellipse is located ad the origin. The smaller ellipses are located
!         on all four sides. The final ellipse is a circle located at the
!         origin.
! =============================================================================
program test_ellipse_multi_merge_1
    use unit_test
    use constants, only : pi, three, four, five
    use parcel_container
    use parcel_merge, only : merge_ellipses
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_ellipse
    implicit none

    double precision :: a1b1, a2b2, error

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call update_parameters()

    call parcel_alloc(5)

    a1b1 = 1.44d0
    a2b2 = 0.25d0


    !
    ! muti-geometric merging
    !

    call parcel_setup

    ! geometric merge
    parcel%lambda = five
    parcel%merge_type = 'multi-geometric'
    parcel%vfraction = three

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse multi-merge 1 (geometric)', error)

    !
    ! muti-optimal merging
    !

    call parcel_setup

    ! optimal merge
    parcel%lambda = five
    parcel%merge_type = 'multi-optimal'
    parcel%vfraction = three

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse multi-merge 1 (optimal)', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            n_parcels = 5
            parcels%position(1, :) = zero
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(1, 2) = zero
            parcels%buoyancy(1) = 1.5d0
            parcels%humidity(1) = 1.3d0

            ! small parcel left
            parcels%position(2, 1) = -1.2d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = a2b2 * pi
            parcels%B(2, 1) = a2b2
            parcels%B(2, 2) = zero
            parcels%buoyancy(2) = 1.8d0
            parcels%humidity(2) = 1.2d0

            ! small parcel right
            parcels%position(3, 1) = 1.2d0
            parcels%position(3, 2) = zero
            parcels%volume(3) = a2b2 * pi
            parcels%B(3, 1) = a2b2
            parcels%B(3, 2) = zero
            parcels%buoyancy(3) = 1.4d0
            parcels%humidity(3) = 1.1d0

            ! small parcel below
            parcels%position(4, 1) = zero
            parcels%position(4, 2) = -1.2d0
            parcels%volume(4) = a2b2 * pi
            parcels%B(4, 1) = a2b2
            parcels%B(4, 2) = zero
            parcels%buoyancy(4) = 1.7d0
            parcels%humidity(4) = 1.0d0

            ! small parcel above
            parcels%position(5, 1) = zero
            parcels%position(5, 2) = 1.2d0
            parcels%volume(5) = a2b2 * pi
            parcels%B(5, 1) = a2b2
            parcels%B(5, 2) = zero
            parcels%buoyancy(5) = 1.5d0
            parcels%humidity(5) = 1.4d0
        end subroutine parcel_setup

        function eval_max_error() result(max_err)
            double precision :: ab, B11, B12, B22, vol, hum, buoy
            double precision :: max_err

            ! reference solution
            ab = a1b1 + four * a2b2  ! a == b since it is a circle
            B11 = ab
            B12 = zero
            B22 = ab
            vol = ab * pi

            buoy = (1.5d0 * a1b1 + (1.8d0 + 1.4d0 + 1.7d0 + 1.5d0) * a2b2) / ab
            hum  = (1.3d0 * a1b1 + (1.2d0 + 1.1d0 + 1.0d0 + 1.4d0) * a2b2) / ab

            max_err = zero

            max_err = max(max_err, abs(dble(n_parcels - 1)))
            max_err = max(max_err, abs(parcels%B(1, 1) - B11))
            max_err = max(max_err, abs(parcels%B(1, 2) - B12))
            max_err = max(max_err, abs(get_B22(parcels%B(1, 1), &
                                               parcels%B(1, 2), &
                                               parcels%volume(1)) - B22))
            max_err = max(max_err, sum(abs(parcels%position(1, :))))
            max_err = max(max_err, abs(parcels%volume(1) - vol))
            max_err = max(max_err, abs(parcels%buoyancy(1) - buoy))
            max_err = max(max_err, abs(parcels%humidity(1) - hum))
        end function eval_max_error



end program test_ellipse_multi_merge_1
