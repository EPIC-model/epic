module taylorgreen
    use parameters, only : flow
    implicit none

    contains
        function get_flow_velocity(pos) result(vel)
            double precision, intent(in) :: pos(2)
            double precision :: xx, yy
            double precision :: vel(2)

            xx = flow%freq(1) * pos(1) + flow%phase(1)
            yy = flow%freq(2) * pos(2) + flow%phase(2)

            vel(1) = flow%amp(1) * cos(xx) * sin(yy)
            vel(2) = flow%amp(2) * sin(xx) * cos(yy)
        end function get_flow_velocity

        ! grad ordering : dudx, dudy, dvdx, dvdy
        function get_flow_gradient(pos) result(grad)
            double precision, intent(in) :: pos(2)
            double precision :: xx, yy
            double precision :: grad(4)

            xx = flow%freq(1) * pos(1) + flow%phase(1)
            yy = flow%freq(2) * pos(2) + flow%phase(2)

            ! du/dx = - a * A * sin(xx) * sin(yy)
            grad(1) = - flow%freq(1) * flow%amp(1) * sin(xx) * sin(yy)

            ! du/dy = b * A * cos(xx) * cos(yy)
            grad(2) = flow%freq(2) * flow%amp(1) * cos(xx) * cos(yy)

            ! dv/dx = a * B * cos(xx) * np.cos(yy)
            grad(3) = flow%freq(1) * flow%amp(2) * cos(xx) * cos(yy)

            ! dv/dy = - b * B * sin(xx) * sin(yy)
            grad(4) = - flow%freq(2) * flow%amp(2) * sin(xx) * sin(yy)

        end function get_flow_gradient

end module taylorgreen
