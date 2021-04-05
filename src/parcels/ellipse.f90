module ellipse
    use constants, only : pi
    implicit none

    contains

        function get_eigenvalue(B11, B12, B22) result(eval)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: eval

            eval = 0.5 * (B11 + B22) + sqrt(0.25 * (B11 - B22) ** 2 + B12 ** 2)
        end function get_eigenvalue

        function get_eigenvector(eval, B12, B22) result(evec)
            double precision, intent(in) :: eval
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: evec(2)

            evec(1) = eval - B22
            evec(2) = B12

            if (abs(evec(1)) + abs(evec(2)) == 0.0d0) then
                evec = evec + epsilon(evec)
            endif

            evec = evec / norm2(evec)

        end function get_eigenvector

        function get_angle(B11, B12) result(angle)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision             :: B22
            double precision             :: eval
            double precision             :: evec(2)
            double precision             :: angle

            B22 = get_B22(B11, B12)

            eval = get_eigenvalue(B11, B12, B22)

            evec = get_eigenvector(eval, B12, B22)

            angle = atan2(evec(2), evec(1))

        end function get_angle

        function get_angles(B, n_parcels) result(angle)
            double precision, intent(in)   :: B(:, :)
            integer,          intent(in)   :: n_parcels
            double precision               :: angle(n_parcels)
            integer                        :: n

            do n = 1, n_parcels
                angle(n) = get_angle(B(n, 1), B(n, 2))
            enddo
        end function get_angles


        elemental function get_B22(B11, B12) result(B22)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision             :: B22

            B22 = (1.0 + B12 ** 2) / B11

        end function get_B22


        function get_ellipse_points(position, volume, B) result(points)
            double precision, intent(in) :: position(2)
            double precision, intent(in) :: volume
            double precision, intent(in) :: B(2)        ! B11, B12
            double precision             :: B22
            double precision             :: c
            double precision             :: eval
            double precision             :: evec(2)
            double precision             :: dx(2)
            double precision             :: points(2, 2)

            B22 = get_B22(B(1), B(2))

            eval = get_eigenvalue(B(1), B(2), B22)

            ! a * b = volume / pi with a and b semi-major and semi-minor axes
            c = sqrt(volume / pi * abs(2.0 * eval - B(1) - B22))

            evec = get_eigenvector(eval, B(2), B22)

            dx = 0.5 * c * evec

            points(1, :) = position - dx
            points(2, :) = position + dx

        end function get_ellipse_points
end module ellipse
