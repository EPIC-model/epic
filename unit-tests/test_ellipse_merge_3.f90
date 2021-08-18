! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test is like test_ellipse_multi_merge_2 but the ellipses
!         are centred at (1.5, 0.2). Hence, this checks periodicity in x.
! =============================================================================
program test_ellipse_multi_merge_3
    use unit_test
    use constants, only : pi, one, two, three, four, five
    use parcel_container
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_ellipse
    use timer
    implicit none

    double precision :: a1b1, a2b2, error

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('parcel merge', merge_timer)

    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call update_parameters

    call parcel_alloc(3)

    a1b1 = 1.44d0
    a2b2 = 0.25d0


    !
    ! muti-geometric merging
    !

    call parcel_setup

    ! geometric merge
    parcel%merge_type = 'geometric'

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error('geometric')

    call print_result_dp('Test ellipse group-merge 3 (geometric)', error)

    !
    ! muti-optimal merging
    !

    call parcel_setup

    ! optimal merge
    parcel%merge_type = 'optimal'

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error('optimal')

    call print_result_dp('Test ellipse group-merge 3 (optimal)', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            double precision :: d

            d = (dsqrt(a1b1) + dsqrt(a2b2)) * f12 * dsqrt(two)

            n_parcels = 3
            parcels%position(1, 1) = 1.5d0
            parcels%position(1, 2) = 0.2d0
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(1, 2) = zero

            ! small parcel left
            parcels%position(2, 1) = 1.5d0 - d
            parcels%position(2, 2) = 0.2d0 - d
            parcels%volume(2) = a2b2 * pi
            parcels%B(2, 1) = a2b2
            parcels%B(2, 2) = zero

            ! small parcel right
            parcels%position(3, 1) = -extent(1) + 1.5d0 + d
            parcels%position(3, 2) = 0.2d0 + d
            parcels%volume(3) = a2b2 * pi
            parcels%B(3, 1) = a2b2
            parcels%B(3, 2) = zero

        end subroutine parcel_setup

        function eval_max_error(method) result(max_err)
            character(*), intent(in) :: method
            double precision :: ab, B11, B12, B22, vol, d, factor, tmp, mu
            double precision :: max_err

            ! reference solution
            d = (dsqrt(a1b1) + dsqrt(a2b2)) * f12 * dsqrt(two)
            ab = a1b1 + two * a2b2
            vol = ab * pi

            B11 = a1b1 ** 2 / ab + two * a2b2 / ab * (four * d ** 2 + a2b2)
            B12 = two * a2b2 / ab * (four * d ** 2)
            B22 = B11

            if (method .eq. 'geometric') then
                factor = ab / dsqrt(B11 * B22 - B12 ** 2)
                B11 = B11 * factor
                B12 = B12 * factor
                B22 = B22 * factor
            else
                tmp = B11
                mu = solve_quartic(tmp, B12, B22, ab)
                B11 = (tmp - mu * B22) / dabs(one - mu ** 2)
                B12 = B12 / dabs(one - mu)
                B22 = (B22 - mu * tmp) / dabs(one - mu ** 2)
            endif

            max_err = zero
            max_err = max(max_err, abs(dble(n_parcels - 1)))
            max_err = max(max_err, abs(parcels%B(1, 1) - B11))
            max_err = max(max_err, abs(parcels%B(1, 2) - B12))
            max_err = max(max_err, abs(get_B22(parcels%B(1, 1), &
                                               parcels%B(1, 2), &
                                               parcels%volume(1)) - B22))
            max_err = max(max_err, sum(abs(parcels%position(1, :) - (/1.5d0, 0.2d0/))))
            max_err = max(max_err, abs(parcels%volume(1) - vol))
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

end program test_ellipse_multi_merge_3
