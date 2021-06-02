! =============================================================================
!                           Module to split ellipses
! =============================================================================
module parcel_split
    use options, only : verbose
    use constants, only : pi, three
    use parameters, only : vcell
    use parcel_container, only : parcel_container_type, n_parcels
    use parcel_ellipse, only : get_eigenvalue      &
                             , get_eigenvector     &
                             , get_B22             &
                             , get_aspect_ratio
    implicit none

    contains

        subroutine split_ellipses(parcels, threshold, vthreshold)
            type(parcel_container_type), intent(inout) :: parcels
            double precision,            intent(in)    :: threshold
            double precision,            intent(in)    :: vthreshold
            double precision                           :: B11
            double precision                           :: B12
            double precision                           :: B22
            double precision                           :: a2, lam, V
            double precision                           :: evec(2)
            double precision                           :: h
            integer                                    :: last_index
            integer                                    :: n

            last_index = n_parcels

            do n = 1, last_index
                B11 = parcels%B(n, 1)
                B12 = parcels%B(n, 2)
                V = parcels%volume(n)
                B22 = get_B22(B11, B12, V)

                a2 = get_eigenvalue(B11, B12, B22)

                ! a/b
                lam = get_aspect_ratio(a2, V)

                if (lam <= threshold .and. V <= vcell/vthreshold) then
                    cycle
                endif

                !
                ! this ellipse is split, i.e., add a new parcel
                !

                evec = get_eigenvector(a2, B11, B12, B22)

                parcels%B(n, 1) = B11 - 0.75d0 * a2 * evec(1) ** 2
                parcels%B(n, 2) = B12 - 0.75d0 * a2 * (evec(1) * evec(2))

                h = 0.25d0 * dsqrt(three * a2)
                parcels%volume(n) = 0.5d0 * V

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1

                parcels%B(n_parcels, :) = parcels%B(n, :)

                parcels%velocity(n_parcels, :) = parcels%velocity(n, :)
                parcels%vorticity(n_parcels, :) = parcels%vorticity(n, :)
                parcels%volume(n_parcels) = parcels%volume(n)
                parcels%buoyancy(n_parcels) = parcels%buoyancy(n)
                parcels%humidity(n_parcels) = parcels%humidity(n)

                parcels%position(n_parcels, :) = parcels%position(n, :) - h * evec
                parcels%position(n, :) = parcels%position(n, :)  + h * evec
            enddo

            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", last_index, "...", n_parcels
            endif

        end subroutine split_ellipses

end module parcel_split
