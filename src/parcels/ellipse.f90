! =============================================================================
!               This module contains all ellipse operations.
! =============================================================================
module ellipse
    use constants, only : pi   &
                        , fpi  &
                        , zero &
                        , two
    implicit none

    contains

        ! the eigenvalue is the semi-major axis squared (a**2)
        function get_eigenvalue(B11, B12, B22) result(a2)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: a2

            a2 = 0.5d0 * (B11 + B22) + dsqrt(0.25d0 * (B11 - B22) ** 2 + B12 ** 2)
        end function get_eigenvalue

        function get_eigenvector(a2, B12, B22) result(evec)
            double precision, intent(in) :: a2
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: evec(2)

            evec(1) = a2 - B22
            evec(2) = B12

            if (abs(evec(1)) + abs(evec(2)) == zero) then
                evec = evec + epsilon(evec)
            endif

            evec = evec / norm2(evec)

        end function get_eigenvector

        function get_angle(B11, B12, volume) result(angle)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: volume
            double precision             :: B22
            double precision             :: a2
            double precision             :: evec(2)
            double precision             :: angle

            B22 = get_B22(B11, B12, volume)

            a2 = get_eigenvalue(B11, B12, B22)

            evec = get_eigenvector(a2, B12, B22)

            angle = atan2(evec(2), evec(1))

        end function get_angle

        function get_angles(B, volume, n_parcels) result(angle)
            double precision, intent(in)   :: B(:, :)
            double precision, intent(in)   :: volume(:, :)
            integer,          intent(in)   :: n_parcels
            double precision               :: angle(n_parcels)
            integer                        :: n

            do n = 1, n_parcels
                angle(n) = get_angle(B(n, 1), B(n, 2), volume(n, 1))
            enddo
        end function get_angles


        elemental function get_B22(B11, B12, volume) result(B22)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: volume
            double precision             :: B22
            double precision             :: ab

            ab = get_ab(volume)

            B22 = (ab ** 2 + B12 ** 2) / B11

        end function get_B22

        elemental function get_ab(volume) result(ab)
            double precision, intent(in) :: volume
            double precision             :: ab

            ab = volume * fpi
        end function

        elemental function get_aspect_ratio(a2, volume) result(lam)
            double precision, intent(in) :: a2
            double precision, intent(in) :: volume
            double precision             :: lam

            lam = (a2 / volume) * pi
        end function

        function get_ellipse_points(position, volume, B) result(points)
            double precision, intent(in) :: position(2)
            double precision, intent(in) :: volume
            double precision, intent(in) :: B(2)        ! B11, B12
            double precision             :: B22
            double precision             :: c
            double precision             :: a2
            double precision             :: evec(2)
            double precision             :: dx(2)
            double precision             :: points(2, 2)

            B22 = get_B22(B(1), B(2), volume)

            a2 = get_eigenvalue(B(1), B(2), B22)

            c = dsqrt(abs(two * a2 - B(1) - B22))

            evec = get_eigenvector(a2, B(2), B22)

            dx = 0.5d0 * c * evec

            points(1, :) = position - dx
            points(2, :) = position + dx

        end function get_ellipse_points
end module ellipse
