module rk4_utils
    use parcel_ellipse, only : get_B22
    use constants, only : two

    implicit none

    contains

        function get_B(Bin, S, volume) result(Bout)
            double precision, intent(in) :: Bin(2)
            double precision, intent(in) :: S(4)
            double precision, intent(in) :: volume
            double precision             :: Bout(2)

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(1) = two * (S(1) * Bin(1) + S(2) * Bin(2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(2) = S(3) * Bin(1) + S(2) * get_B22(Bin(1), Bin(2), volume)

        end function get_B

        function get_stretch(S, n) result(l)
            use parcel_point, only : get_eigenvalue
            double precision, intent(in) :: S(:, :)
            integer,          intent(in) :: n
            integer                      :: i
            double precision             :: l(n)

            do i = 1, n
                l(i) = get_eigenvalue(S(i, :))
            enddo
        end function get_stretch

end module rk4_utils
