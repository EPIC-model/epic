module straka
    use physics
    use options, only : parcel_info
    use constants
    use parcel_container, only : parcels, n_parcels
    use fields, only : get_position
    implicit none

    contains

        subroutine straka_init
            integer          :: n
            double precision :: xc, xr, zc, zr, L, dT, dtheta, theta_0

            ! in metres
            xc = zero
            xr = four
            zc = three
            zr = two

            ! reference potential temperature
            theta_0 = 300.0d0 ! Kelvin

            dT = zero
            do n = 1, n_parcels
                L = ((parcels%position(n, 1) - xc) / xr) ** 2 &
                  + ((parcels%position(n, 2) - zc) / zr) ** 2

                L = dsqrt(L)

                ! temperature perturbation
                if (L > one) then
                    dT = zero
                else
                    dT = -7.5d0 * (dcos(pi * L) + one)
                endif

                ! potential temperatur
                dtheta = dT * pi

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                ! (dtheta = theta - theta_0)
                parcels%buoyancy(n) = gravity * dtheta / theta_0
            enddo
        end subroutine straka_init
end module
