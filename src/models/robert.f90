! =============================================================================
!                       Robert ("Rising bubble") test case
!
!   This module sets up a bubble with uniform or Gaussian buoyancy profile.
!   The setup is from the paper
!
!   Robert, A. (1993). Bubble Convection Experiments with a Semi-implicit
!   Formulation of the Euler Equations, Journal of Atmospheric Sciences,
!   50(13), 1865-1873. Retrieved Jun 2, 2021, from
!   https://doi.org/10.1175/1520-0469(1993)050<1865:BCEWAS>2.0.CO;2
!
! =============================================================================

module robert
    use physics
    use options, only : parcel_info, robert_opt
    use constants
    use parcel_container, only : parcels, n_parcels
    use fields, only : get_position
    implicit none

    private :: robert_uniform_init, &
               robert_gaussian_init

    contains

        subroutine robert_init
            character(:), allocatable :: distr

            distr = trim(robert_opt%distr)

            select case (distr)
                case ('uniform')
                    call robert_uniform_init
                case ('gaussian')
                    call robert_gaussian_init
                case default
                    print *, "Invalid distribution type: '", distr, "'"
                    stop
            end select

        end subroutine robert_init

        subroutine robert_uniform_init
            integer          :: n
            double precision :: xc, zc, r, r2, dtheta, theta_ref, pos(2)

            ! in metres
            xc = robert_opt%center(1)
            zc = robert_opt%center(2)
            r2 = robert_opt%outer_radius ** 2

            ! reference potential temperature
            theta_ref = robert_opt%theta_ref

            do n = 1, n_parcels
                r = (parcels%position(n, 1) - xc) ** 2 &
                  + (parcels%position(n, 2) - zc) ** 2

                ! potential temperature perturbation
                dtheta = zero

                if (r <= r2) then
                    dtheta = robert_opt%theta_max
                endif

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                ! (dtheta = theta - theta_ref)
                parcels%buoyancy(n) = gravity * dtheta / theta_ref
            enddo
        end subroutine robert_uniform_init


        subroutine robert_gaussian_init
            integer          :: n
            double precision :: xc, zc, r, rr, dtheta, theta_ref, theta_max, pos(2), a, s, fs2

            ! in metres
            xc = robert_opt%center(1)
            zc = robert_opt%center(2)
            rr = robert_opt%outer_radius
            a  = robert_opt%inner_radius
            s  = robert_opt%width
            fs2 = (one / s) ** 2

            ! reference potential temperature
            theta_ref = robert_opt%theta_ref

            ! max. potential temperature perturbation
            theta_max = robert_opt%theta_max

            do n = 1, n_parcels
                r = (parcels%position(n, 1) - xc) ** 2 &
                  + (parcels%position(n, 2) - zc) ** 2

                r = dsqrt(r)

                ! potential temperature perturbation
                dtheta = zero

                if (r <= rr) then
                    if (r <= a) then
                        dtheta = theta_max
                    else
                        dtheta = theta_max * dexp(-(r - a) ** 2 * fs2)
                    endif
                endif

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                ! (dtheta = theta - theta_ref)
                parcels%buoyancy(n) = gravity * dtheta / theta_ref
            enddo
        end subroutine robert_gaussian_init
end module
