! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the merging of three ellipses. The biggest
!         ellipse is located ad the origin. The smaller ellipses are located
!         tangentially at 45 and 225 degrees. The final ellipse is an ellipse
!         located at the origin.
! =============================================================================
program test_ellipse_multi_merge_2
    use unit_test
    use constants, only : pi, one, two, four
    use parcel_container
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, nx, nz, lower, extent
    use parcel_ellipse
    use timer
    implicit none

    double precision :: a1b1, a2b2, error

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('parcel merge', merge_timer)

    call update_parameters

    call parcel_alloc(3)

    a1b1 = 1.44d0
    a2b2 = 0.25d0


    !
    ! muti-geometric merging
    !

    call parcel_setup

    ! geometric merge
    parcel%lambda = 5.0
    parcel%merge_type = 'multi-geometric'
    parcel%vfraction = 3

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error('multi-geometric')

    call print_result_dp('Test ellipse multi-merge 2 (geometric)', error)

    !
    ! muti-optimal merging
    !

    call parcel_setup

    ! optimal merge
    parcel%lambda = 5.0
    parcel%merge_type = 'multi-optimal'
    parcel%vfraction = 3

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error('multi-optimal')

    call print_result_dp('Test ellipse multi-merge 2 (optimal)', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            double precision :: d

            d = (dsqrt(a1b1) + dsqrt(a2b2)) * f12 * dsqrt(two)

            n_parcels = 3
            parcels%position(1, :) = zero
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(1, 2) = zero
            parcels%buoyancy(1) = 1.5d0
            parcels%humidity(1) = 1.3d0

            ! small parcel left
            parcels%position(2, 1) = -d
            parcels%position(2, 2) = -d
            parcels%volume(2) = a2b2 * pi
            parcels%B(2, 1) = a2b2
            parcels%B(2, 2) = zero
            parcels%buoyancy(2) = 1.8d0
            parcels%humidity(2) = 1.2d0

            ! small parcel right
            parcels%position(3, 1) = d
            parcels%position(3, 2) = d
            parcels%volume(3) = a2b2 * pi
            parcels%B(3, 1) = a2b2
            parcels%B(3, 2) = zero
            parcels%buoyancy(3) = 1.4d0
            parcels%humidity(3) = 1.1d0

        end subroutine parcel_setup

        function eval_max_error(method) result(max_err)
            character(*), intent(in) :: method
            double precision :: ab, B11, B12, B22, vol, angle, d, factor, tmp, mu
            double precision :: max_err, hum, buoy

            ! reference solution
            d = (dsqrt(a1b1) + dsqrt(a2b2)) * f12 * dsqrt(two)
            ab = a1b1 + two * a2b2
            vol = ab * pi

            buoy = (1.5d0 * a1b1 + (1.8d0 + 1.4d0) * a2b2) / ab
            hum  = (1.3d0 * a1b1 + (1.2d0 + 1.1d0) * a2b2) / ab

            B11 = a1b1 ** 2 / ab + two * a2b2 / ab * (four * d ** 2 + a2b2)
            B12 = two * a2b2 / ab * (four * d ** 2)
            B22 = B11

            if (method .eq. 'multi-geometric') then
                factor = ab / dsqrt(B11 * B22 - B12 ** 2)
                B11 = B11 * factor
                B12 = B12 * factor
                B22 = B22 * factor
            else
                tmp = B11
                mu = solve_quartic(tmp, B12, B22, ab)
                B11 = (tmp - mu * B22) / (one - mu ** 2)
                B12 = B12 / (one - mu)
                B22 = (B22 - mu * tmp) / (one - mu ** 2)
            endif

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


        function solve_quartic(B11, B12, B22, ab) result(mu)
            double precision, intent(in) :: B11, B12, B22, ab
            double precision             :: mu, detB, merr, mup
            double precision             :: a, b ,c, a2b2i
            ! Solve the quartic to find best fit ellipse:
            !
            !
            !   Newton-Raphson to get smallest root:
            !      mu_{n+1} = mu_{n} - f(mu_{n}) / f'(mu_{n})
            !
            !      where
            !          f(mu_{n})  = mu_{n} ** 4 + b * mu_{n} ** 2 + a * mu_{n} + c
            !          f'(mu_{n}) = 4 * mu_{n} ** 3 + b * 2 * mu_{n} + a
            a2b2i = one / ab ** 2
            detB = (B11 * B22 - B12 * B12) * a2b2i
            a = (B11 ** 2 + B22 ** 2 + two * B12 ** 2) * a2b2i
            b = -two - detB
            c = one - detB

            ! initial guess
            mu = - c / a

            merr = 1.0
            do while (merr > 1.e-12)
                mup = (c + mu * (a + mu * (b + mu * mu))) / (a + mu * (two * b + four * mu * mu))
                mu = mu - mup
                merr = abs(mup)
            enddo
        end function solve_quartic

end program test_ellipse_multi_merge_2
