module taylorgreen
    use parameters, only : amp, freq, phase
    implicit none

    contains
        function get_flow_velocity(x, y) result(vel)
            double precision, intent(in) :: x, y
            double precision :: xx, yy
            double precision :: vel(2)

            xx = amp(1) * x + phase(1)
            yy = amp(2) * y + phase(2)

            vel(1) = amp(1) * cos(xx) * sin(yy)
            vel(2) = amp(2) * sin(xx) * cos(yy)
        end function get_flow_velocity

        ! grad ordering : dudx, dudy, dvdx, dvdy
        function get_flow_gradient(x, y) result(grad)
            double precision, intent(in) :: x, y
            double precision :: xx, yy
            double precision :: grad(4)

            xx = amp(1) * x + phase(1)
            yy = amp(2) * y + phase(2)

            ! du/dx = - a * A * sin(xx) * sin(yy)
            grad(1) = - freq(1) * amp(1) * sin(xx) * sin(yy)

            ! du/dy = b * A * cos(xx) * cos(yy)
            grad(2) = freq(2) * amp(1) * cos(xx) * cos(yy)

            ! dv/dx = a * B * cos(xx) * np.cos(yy)
            grad(3) = freq(1) * amp(2) * cos(xx) * cos(yy)

            ! dv/dy = - b * B * sin(xx) * sin(yy)
            grad(4) = - freq(2) * amp(2) * sin(xx) * sin(yy)

        end function get_flow_gradient

end module taylorgreen
