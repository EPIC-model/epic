module ellipse
    use constants, only : pi
    use parcel_container, only : parcel_container_type, n_parcels
    implicit none

    private :: get_eigenvalue,  &
               get_eigenvector, &
               get_B22

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

            ! add epsilon if eval == B22 --> circle
            if (abs(evec(1)) <= epsilon(eval)) then
                evec(1) = evec(1) + epsilon(evec(1))
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


        function get_B22(B11, B12) result(B22)
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision             :: B22

            B22 = (1.0 + B12 ** 2) / B11

        end function get_B22


        function get_ellipse_points(position, volume, B11, B12) result(points)
            double precision, intent(in) :: position
            double precision, intent(in) :: volume
            double precision, intent(in) :: B11
            double precision, intent(in) :: B12
            double precision             :: B22
            double precision             :: c
            double precision             :: eval
            double precision             :: evec(2)
            double precision             :: dx(2)
            double precision             :: points(2, 2)

            B22 = get_B22(B11, B12)

            eval = get_eigenvalue(B11, B12, B22)

            ! a * b = volume / pi with a and b semi-major and semi-minor axes
            c = sqrt(volume / pi * abs(2.0 * eval - B11 - B22))

            evec = get_eigenvector(eval, B12, B22)

            dx = 0.5 * c * evec

            points(:, 1) = position - dx
            points(:, 2) = position + dx

        end function get_ellipse_points

        subroutine split_ellipse(parcels, threshold)
            type(parcel_container_type), intent(inout) :: parcels
            double precision,            intent(in)    :: threshold
            double precision                           :: B11
            double precision                           :: B12
            double precision                           :: B22
            double precision                           :: eval
            double precision                           :: evec(2)
            double precision                           :: h
            integer                                    :: last_index
            integer                                    :: i

            last_index = n_parcels

            do i = 1, last_index
                B11 = parcels%B11(i)
                B12 = parcels%B12(i)
                B22 = get_B22(B11, B12)

                eval = get_eigenvalue(B11, B12, B22)

                if (eval > threshold) then
                    cycle
                endif

                !
                ! this ellipse is split, i.e., add a new parcel
                !

                evec = get_eigenvector(eval, B12, B22)

                parcels%B11(i) = 2.0 * B11 - 1.5 * eval * evec(1) ** 2
                parcels%B12(i) = 2.0 * B12 - 1.5 * eval * (evec(1) * evec(2))

                h = 0.25 * sqrt(3.0 * eval * parcels%volume(i) / pi)
                parcels%volume(i) = 0.5 * parcels%volume(i)

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1

                parcels%velocity(n_parcels, :) = parcels%velocity(i, :)
                parcels%volume(n_parcels) = parcels%volume(i)

                parcels%position(n_parcels, :) = parcels%position(i, :) - h * evec
                parcels%position(i, :) = parcels%position(i, :)  + h * evec
            enddo
        end subroutine split_ellipse

end module ellipse
