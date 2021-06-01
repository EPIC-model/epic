module robert
    use physics
    use options, only : parcel_info
    use constants
    use parcel_container, only : parcels, n_parcels
    use fields, only : get_position
    implicit none

    contains

        subroutine robert_uniform_init
            integer          :: n
            double precision :: xc, zc, r, r2, dT, dtheta, theta_0, pos(2)

            ! in metres
            xc = zero
            zc = 250.0d0 + ten
            r2 = 62500.d0

            ! reference potential temperature
            theta_0 = 303.15d0 * pi ! Kelvin (approx 30 degree Celsius)

            dT = zero
            do n = 1, n_parcels
                r = (parcels%position(n, 1) - xc) ** 2 &
                  + (parcels%position(n, 2) - zc) ** 2

                ! temperature perturbation
                dT = zero

                if (r <= r2) then
                    dT = 0.5d0
                endif

                ! potential temperatur
                dtheta = dT * pi

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                ! (dtheta = theta - theta_0)
                parcels%buoyancy(n) = gravity * dtheta / theta_0
            enddo
        end subroutine robert_uniform_init


        subroutine robert_gaussian_init
            integer          :: n
            double precision :: xc, zc, r, rr, dT, dtheta, theta_0, pos(2), a, fs2

            ! in metres
            xc  = zero
            zc  = 250.0d0 + ten
            rr  = 250.0d0
            a   = 50.0d0
            fs2 = (one / hundred) ** 2

            ! reference potential temperature
            theta_0 = 303.15d0 * pi ! Kelvin (approx 30 degree Celsius)

            dT = zero
            do n = 1, n_parcels
                r = (parcels%position(n, 1) - xc) ** 2 &
                  + (parcels%position(n, 2) - zc) ** 2

                r = dsqrt(r)

                ! temperature perturbation
                dT = zero

                if (r <= rr) then
                    if (r <= a) then
                        dT = 0.5d0
                    else
                        dT = 0.5d0 * dexp(-(r - a) ** 2 * fs2)
                    endif
                endif

                ! potential temperatur
                dtheta = dT * pi

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                ! (dtheta = theta - theta_0)
                parcels%buoyancy(n) = gravity * dtheta / theta_0
            enddo
        end subroutine robert_gaussian_init


end module
