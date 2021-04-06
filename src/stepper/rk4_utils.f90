module rk4_utils
    use constants, only : max_num_parcels
    use ellipse, only : get_B22
    use parcel_container

    implicit none

    contains

        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(max_num_parcels, 2)
            double precision, intent(in) :: S(max_num_parcels, 4)
            double precision             :: Bout(max_num_parcels, 2)

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(1:n_parcels, 1) = 2.0 * (S(1:n_parcels, 1) * Bin(1:n_parcels, 1) &
                                 + S(1:n_parcels, 2) * Bin(1:n_parcels, 2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(1:n_parcels, 2) = S(1:n_parcels, 3) * Bin(1:n_parcels, 1) &
                                 + S(1:n_parcels, 2) * get_B22(Bin(1:n_parcels, 1), Bin(1:n_parcels, 2))

        end function get_B

end module rk4_utils
