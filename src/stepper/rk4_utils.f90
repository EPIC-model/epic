module rk4_utils
    use ellipse, only : get_B22
    use parcel_container

    implicit none

    contains

        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(:, :)
            double precision, intent(in) :: S(:, :)
            double precision, dimension(size(Bin,1),size(Bin,2)) :: Bout

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(:, 1) = 2.0 * (S(:, 1) * Bin(:, 1) &
                                 + S(:, 2) * Bin(:, 2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(:, 2) = S(:, 3) * Bin(:, 1) &
                                 + S(:, 2) * get_B22(Bin(:, 1), Bin(:, 2))

        end function get_B

end module rk4_utils
