! =============================================================================
!                       Straka ("Cold bubble") test case
!
!   This module sets up a cold bubble as described in Straka et al (1993).
!
!   Straka, J.M., Wilhelmson, R.B., Wicker, L.J., Anderson, J.R. and
!   Droegemeier, K.K. (1993), Numerical solutions of a non-linear density
!   current: A benchmark solution and comparisons. Int. J. Numer. Meth. Fluids,
!   17: 1-22. https://doi.org/10.1002/fld.1650170103
!
! =============================================================================

module straka
    use physics
    use options, only : parcel_info, straka_opt
    use constants
    use parcel_container, only : parcels, n_parcels
    use fields, only : get_position
    implicit none

    contains

        subroutine straka_init
            integer          :: n
            double precision :: xc, xr, zc, zr, L, dtheta, theta_ref, theta_max, pos(2)

            ! in metres
            xc = straka_opt%center(1)
            xr = straka_opt%radii(2)
            zc = straka_opt%center(2)
            zr = straka_opt%radii(2)

            ![Kelvin] reference potential temperature
            theta_ref = straka_opt%theta_ref

            theta_max = -0.5d0 * straka_opt%theta_max

            do n = 1, n_parcels
                L = ((parcels%position(n, 1) - xc) / xr) ** 2 &
                  + ((parcels%position(n, 2) - zc) / zr) ** 2

                L = dsqrt(L)

                ! potential temperature perturbation
                dtheta = zero

                if (L <= one) then
                    dtheta = theta_max * (dcos(pi * L) + one)
                endif

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                ! (dtheta = theta - theta_ref)
                parcels%buoyancy(n) = gravity * dtheta / theta_ref
            enddo
        end subroutine straka_init
end module
