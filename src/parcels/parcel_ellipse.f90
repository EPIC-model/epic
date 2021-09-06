! =============================================================================
!               This module contains all ellipse operations.
! =============================================================================
module parcel_ellipse
    use constants, only : pi   &
                        , fpi  &
                        , zero &
                        , two  &
                        , f12  &
                        , f14
    implicit none

    contains

        ! Obtain the largest eigenvalue (i.e. the semi-major axis squared [a**2])
        ! @param[in] B11 is the B matrix (1, 1) element
        ! @param[in] B12 is the B matrix (1, 2) element
        ! @param[in] B22 is the B matrix (2, 2) element
        ! @returns the largest eigenvalue
        function get_eigenvalue(B11, B12, B22) result(a2)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: a2

            a2 = f12 * (B11 + B22) + dsqrt(f14 * (B11 - B22) ** 2 + B12 ** 2)
        end function get_eigenvalue

        ! Obtain the eigenvector of the largest eigenvalue
        ! @param[in] a2 is the largest eigenvalue
        ! @param[in] B11 is the B matrix (1, 1) element
        ! @param[in] B12 is the B matrix (1, 2) element
        ! @param[in] B22 is the B matrix (2, 2) element
        ! @returns the eigenvector
        function get_eigenvector(a2, B11, B12, B22) result(evec)
            double precision, intent(in) :: a2
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: B22
            double precision             :: evec(2)

            evec(1) = a2 - B22
            evec(2) = B12

            if (dabs(evec(1)) + dabs(evec(2)) == zero) then
                if (B11 > B22) then
                    evec(1) = evec(1) + epsilon(evec(1))
                else
                    evec(2) = evec(2) + epsilon(evec(2))
                endif
            endif

            evec = evec / norm2(evec)

        end function get_eigenvector

        ! used in unit tests only
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

            evec = get_eigenvector(a2, B11, B12, B22)

            angle = datan2(evec(2), evec(1))

        end function get_angle

        ! Obtain the B22 matrix element
        ! @param[in] B11 is the B matrix (1, 1) element
        ! @param[in] B12 is the B matrix (1, 2) element
        ! @param[in] volume of the parcel
        ! @returns B22
        elemental function get_B22(B11, B12, volume) result(B22)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision, intent(in) :: volume
            double precision             :: B22
            double precision             :: ab

            ab = get_ab(volume)

            B22 = (ab ** 2 + B12 ** 2) / B11

        end function get_B22

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] volume of the parcel
        ! @returns ab = volume / pi
        elemental function get_ab(volume) result(ab)
            double precision, intent(in) :: volume
            double precision             :: ab

            ab = volume * fpi
        end function

        ! Obtain the aspect ratio of the parcel(s).
        ! @param[in] a2 is the largest eigenvalue
        ! @param[in] volume of the parcel(s)
        ! @returns lam = a/b
        elemental function get_aspect_ratio(a2, volume) result(lam)
            double precision, intent(in) :: a2
            double precision, intent(in) :: volume
            double precision             :: lam

            lam = (a2 / volume) * pi
        end function

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] position vector of the parcel
        ! @param[in] volume of the parcel
        ! @param[in] B matrix elements of the parcel
        ! @returns the parcel support points
        function get_ellipse_points(position, volume, B) result(points)
            double precision, intent(in) :: position(2)
            double precision, intent(in) :: volume
            double precision, intent(in) :: B(2)        ! B11, B12
            double precision             :: B22
            double precision             :: c
            double precision             :: a2
            double precision             :: evec(2)
            double precision             :: h(2)
            double precision             :: points(2, 2)

            B22 = get_B22(B(1), B(2), volume)

            a2 = get_eigenvalue(B(1), B(2), B22)

            c = dsqrt(dabs(two * a2 - B(1) - B22))

            evec = get_eigenvector(a2, B(1), B(2), B22)

            h = f12 * c * evec

            points(1, :) = position - h
            points(2, :) = position + h

        end function get_ellipse_points
end module parcel_ellipse
