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
        function get_eigenvalue(B) result(a2)
            double precision, intent(in) :: B(3)
            double precision             :: a2

            a2 = f12 * (B(1) + B(3)) + dsqrt(f14 * (B(1) - B(3)) ** 2 + B(2) ** 2)
        end function get_eigenvalue

        ! Obtain the eigenvector of the largest eigenvalue
        ! @param[in] a2 is the largest eigenvalue
        ! @param[in] B11 is the B matrix (1, 1) element
        ! @param[in] B12 is the B matrix (1, 2) element
        ! @param[in] B22 is the B matrix (2, 2) element
        ! @returns the eigenvector
        function get_eigenvector(a2, B) result(evec)
            double precision, intent(in) :: a2
            double precision, intent(in) :: B(3)
            double precision             :: evec(2)

            evec(1) = a2 - B(3)
            evec(2) = B(2)

            if (dabs(evec(1)) + dabs(evec(2)) == zero) then
                if (B(1) > B(3)) then
                    evec(1) = evec(1) + epsilon(evec(1))
                else
                    evec(2) = evec(2) + epsilon(evec(2))
                endif
            endif

            evec = evec / norm2(evec)

        end function get_eigenvector

        ! used in unit tests only
        function get_angle(B) result(angle)
            double precision, intent(in) :: B(3)
            double precision             :: a2
            double precision             :: evec(2)
            double precision             :: angle

            a2 = get_eigenvalue(B)

            evec = get_eigenvector(a2, B)

            angle = datan2(evec(2), evec(1))

        end function get_angle

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] area of the parcel
        ! @returns ab = area / pi
        elemental function get_ab(area) result(ab)
            double precision, intent(in) :: area
            double precision             :: ab

            ab = area * fpi
        end function

        function get_area(B) result(area)
            double precision, intent(in) :: B(3) ! = (B11, B12, B22)
            double precision             :: area

#ifndef NDEBUG
            if (B(1) * B(3) < B(2) ** 2) then
                print *, "Error: Determinant is negative. Unable to calculate parcel area."
                stop
            endif
#endif
            area = pi * dsqrt(B(1) * B(3) - B(2) ** 2)
        end function get_area

        ! Obtain the aspect ratio a/b of the parcel(s).
        ! @param[in] a2 is the largest eigenvalue
        ! @param[in] area of the parcel(s)
        ! @returns lam = a/b
        elemental function get_aspect_ratio(a2, area) result(lam)
            double precision, intent(in) :: a2
            double precision, intent(in) :: area
            double precision             :: lam

            lam = (a2 / area) * pi
        end function

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] position vector of the parcel
        ! @param[in] B matrix elements of the parcel
        ! @returns the parcel support points
        function get_ellipse_points(position, B) result(points)
            double precision, intent(in) :: position(2)
            double precision, intent(in) :: B(3)        ! B11, B12, B22
            double precision             :: c
            double precision             :: a2
            double precision             :: evec(2)
            double precision             :: h(2)
            double precision             :: points(2, 2)

            a2 = get_eigenvalue(B)

            c = dsqrt(dabs(two * a2 - B(1) - B(3)))

            evec = get_eigenvector(a2, B)

            h = f12 * c * evec

            points(:, 1) = position - h
            points(:, 2) = position + h

        end function get_ellipse_points
end module parcel_ellipse
