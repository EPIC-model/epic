module parcel_merge
    use constants, only : pi
    use parcel_container, only : parcel_container_type
    use ellipse, only : get_B22
    implicit none

    private :: geometric_merge

    contains
        subroutine merge_ellipses(parcels)
            type(parcel_container_type), intent(inout) :: parcels

            ! find parcels to merge


        end subroutine merge_ellipses

        ! merge ith parcel into jth parcel
        subroutine geometric_merge(parcels, n, m)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: n ! index of smaller parcel
            integer,                     intent(in)    :: m ! index of larger parcel
            double precision                           :: B11_1, B11_2
            double precision                           :: B12_1, B12_2
            double precision                           :: B22_1, B22_2, B22
            double precision                           :: a1b1, a2b2, ab, isqrab
            double precision                           :: mu1, mu2, zet, eta
            double precision                           :: mu11, mu22, mu12, detB

            B11_1 = parcels%B(n, 1)
            B11_2 = parcels%B(m, 1)

            B12_1 = parcels%B(n, 2)
            B12_2 = parcels%B(m, 2)

            B22_1 = get_B22(B11_1, B12_1)
            B22_2 = get_B22(B11_2, B12_2)

            a1b1 = parcels%volume(n, 1) / pi
            a2b2 = parcels%volume(m, 1) / pi

            ab = a1b1 + a2b2
            isqrab = 1.0 / sqrt(ab)

            mu1 = a1b1 / ab
            mu2 = a2b2 / ab

            zet = 2.0 * isqrab * (parcels%position(m, 1) - parcels%position(n, 1))
            eta = 2.0 * isqrab * (parcels%position(m, 2) - parcels%position(n, 2))

            mu12 = mu1 * mu2
            mu11 = mu1 * mu1
            mu22 = mu2 * mu2

            parcels%B(m, 1) = mu12 * zet ** 2  + mu11 * B11_1 + mu22 * B11_2
            parcels%B(m, 2) = mu12 * zet * eta + mu11 * B12_1 + mu22 * B12_2
            B22             = mu12 * eta ** 2  + mu11 * B22_2 + mu22 * B22_2

            ! normalize such that determinant of the merger is 1
            detB = parcels%B(m, 1) * B22 - parcels%B(m, 2) ** 2

            parcels%B(m, 1) = parcels%B(m, 1) / detB
            parcels%B(m, 2) = parcels%B(m, 2) / detB

        end subroutine geometric_merge


        subroutine find_nearest

        end subroutine find_nearest



end module parcel_merge
