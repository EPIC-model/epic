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
    use options, only : parcel_info, robert_opt, bubble_type
    use constants
    use parcel_container, only : parcels, n_parcels
    use fields, only : get_position
    implicit none

    private :: robert_uniform_init, &
               robert_gaussian_init

    contains

        subroutine robert_init
            integer           :: k
            type(bubble_type) :: bubble

            parcels%buoyancy(1:n_parcels) = zero

            if (robert_opt%n_bubbles > size(robert_opt%bubbles)) then
                print *, 'Number of bubbles beyond upper limit.'
                stop
            endif


            do k = 1, robert_opt%n_bubbles
                bubble = robert_opt%bubbles(k)

                select case (trim(bubble%distr))
                    case ('uniform')
                        call robert_uniform_init(bubble)
                    case ('gaussian')
                        call robert_gaussian_init(bubble)
                    case default
                        print *, "Invalid distribution type: '", trim(bubble%distr), "'"
                        stop
                end select
            enddo
        end subroutine robert_init

        subroutine robert_uniform_init(bubble)
            type(bubble_type), intent(in) :: bubble
            integer                       :: n
            double precision              :: xc, zc, r2, r2_out, dtheta, theta_ref, pos(2)

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r2_out = bubble%outer_radius ** 2

            ! reference potential temperature
            theta_ref = robert_opt%theta_ref

            do n = 1, n_parcels
                r2 = (parcels%position(n, 1) - xc) ** 2 &
                   + (parcels%position(n, 2) - zc) ** 2

                ! potential temperature perturbation
                dtheta = zero

                if (r2 <= r2_out) then
                    dtheta = bubble%theta_max
                endif

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                ! (dtheta = theta - theta_ref)
                parcels%buoyancy(n) = parcels%buoyancy(n) &
                                    + gravity * dtheta / theta_ref
            enddo
        end subroutine robert_uniform_init


        subroutine robert_gaussian_init(bubble)
            type(bubble_type), intent(in) :: bubble
            integer                       :: n
            double precision              :: xc, zc, r, r_out, r_in, dtheta
            double precision              :: theta_ref, theta_max, pos(2), s, fs2

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r_out = bubble%outer_radius
            r_in = bubble%inner_radius
            s = bubble%width
            fs2 = (one / s) ** 2

            ! reference potential temperature
            theta_ref = robert_opt%theta_ref

            ! max. potential temperature perturbation
            theta_max = bubble%theta_max

            do n = 1, n_parcels
                r = (parcels%position(n, 1) - xc) ** 2 &
                  + (parcels%position(n, 2) - zc) ** 2

                r = dsqrt(r)

                ! potential temperature perturbation
                dtheta = zero

                if (r <= r_out) then
                    if (r <= r_in) then
                        dtheta = theta_max
                    else
                        dtheta = theta_max * dexp(-(r - r_in) ** 2 * fs2)
                    endif
                endif

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                ! (dtheta = theta - theta_ref)
                parcels%buoyancy(n) = parcels%buoyancy(n) &
                                    + gravity * dtheta / theta_ref
            enddo
        end subroutine robert_gaussian_init
end module
