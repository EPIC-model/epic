module rk4_utils
    use constants, only : max_num_parcels
    use ellipse, only : get_B22

    implicit none

    contains

        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(max_num_parcels, 2)
            double precision, intent(in) :: S(max_num_parcels, 4)
            double precision             :: Bout(max_num_parcels, 2)
            double precision             :: B22(max_num_parcels)

            B22 = get_B22(Bin(:, 1), Bin(:, 2))

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(:, 1) = 2.0 * (S(:, 1) * Bin(:, 1) + S(:, 2) * Bin(:, 2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(:, 2) = S(:, 3) * Bin(:, 1) + S(:, 2) * B22

        end function get_B

end module rk4_utils
