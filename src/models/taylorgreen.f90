! =============================================================================
!                           2D TAYLOR-GREEN FLOW
!
!                   u(x, y) = A * cos(ax + d) * sin(by + e)
!                   v(x, y) = B * sin(ax + d) * cos(by + e)
!
!                   The code uses following variable names:
!                       taylor_green_opt%freq  = (/a, b/)
!                       taylor_green_opt%amp   = (/A, B/)
!                       taylor_green_opt%phase = (/d, e/)
! =============================================================================
module taylorgreen
    use options, only : taylor_green_opt
    use constants, only : pi
    use parcel_container, only : parcels, n_parcels
    use fields
    implicit none

    private
        double precision, parameter :: hpi = 0.5d0 * pi

    public :: get_flow_velocity,    &
              get_flow_gradient,    &
              get_flow_vorticity,   &
              taylorgreen_init

    contains
        subroutine taylorgreen_init
            integer          :: n

            do n = 1, n_parcels
                parcels%vorticity(n, :) = get_flow_vorticity(parcels%position(n, :))
            enddo

        end subroutine taylorgreen_init


        function get_flow_velocity(pos) result(vel)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: vel(2)

            call get_flow_pos(pos, xx, zz)

            vel(1) = taylor_green_opt%amp(1) * dcos(xx) * dsin(zz)
            vel(2) = taylor_green_opt%amp(2) * dsin(xx) * dcos(zz)
        end function get_flow_velocity

        ! grad ordering : dudx, dudy, dvdx, dvdy
        function get_flow_gradient(pos) result(grad)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: grad(4)

            call get_flow_pos(pos, xx, zz)

            ! du/dx = - a * A * sin(xx) * sin(zz)
            grad(1) = - taylor_green_opt%freq(1) * taylor_green_opt%amp(1) * dsin(xx) * dsin(zz)

            ! du/dy = b * A * cos(xx) * cos(zz)
            grad(2) = taylor_green_opt%freq(2) * taylor_green_opt%amp(1) * dcos(xx) * dcos(zz)

            ! dv/dx = a * B * cos(xx) * np.cos(zz)
            grad(3) = taylor_green_opt%freq(1) * taylor_green_opt%amp(2) * dcos(xx) * dcos(zz)

            ! dv/dy = - b * B * sin(xx) * sin(zz)
            grad(4) = - taylor_green_opt%freq(2) * taylor_green_opt%amp(2) * dsin(xx) * dsin(zz)

        end function get_flow_gradient

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: omega

            call get_flow_pos(pos, xx, zz)

            omega = (taylor_green_opt%amp(2) * taylor_green_opt%freq(1)     &
                   - taylor_green_opt%amp(1) * taylor_green_opt%freq(2))    &
                   * dcos(xx) * dcos(zz)
        end function get_flow_vorticity

        subroutine get_flow_pos(pos, xx, zz)
            double precision, intent(in) :: pos(2)
            double precision, intent(out) :: xx, zz

            xx = taylor_green_opt%freq(1) * pos(1) + taylor_green_opt%phase(1) + hpi
            zz = taylor_green_opt%freq(2) * pos(2) + taylor_green_opt%phase(2)

        end subroutine get_flow_pos

end module taylorgreen
