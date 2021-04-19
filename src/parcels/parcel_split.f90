! =============================================================================
!                           Module to split ellipses
! =============================================================================
module parcel_split
    use options, only : verbose
    use constants, only : pi
    use parcel_container, only : parcel_container_type, n_parcels
    use ellipse, only : get_eigenvalue, get_eigenvector, get_B22
    implicit none

    contains

        subroutine split_ellipses(parcels, threshold)
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
                B11 = parcels%B(i, 1)
                B12 = parcels%B(i, 2)
                B22 = get_B22(B11, B12)

                eval = get_eigenvalue(B11, B12, B22)

                if (eval <= threshold) then
                    cycle
                endif

                !
                ! this ellipse is split, i.e., add a new parcel
                !

                evec = get_eigenvector(eval, B12, B22)

                parcels%B(i, 1) = 2.0 * B11 - 1.5 * eval * evec(1) ** 2
                parcels%B(i, 2) = 2.0 * B12 - 1.5 * eval * (evec(1) * evec(2))

                h = 0.25 * sqrt(3.0 * eval * parcels%volume(i, 1) / pi)
                parcels%volume(i, 1) = 0.5 * parcels%volume(i, 1)

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1

                parcels%B(n_parcels, :) = parcels%B(i, :)

                parcels%velocity(n_parcels, :) = parcels%velocity(i, :)
                parcels%volume(n_parcels, 1) = parcels%volume(i, 1)

                parcels%position(n_parcels, :) = parcels%position(i, :) - h * evec
                parcels%position(i, :) = parcels%position(i, :)  + h * evec
            enddo

            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", last_index, "...", n_parcels
            endif

        end subroutine split_ellipses

end module parcel_split
