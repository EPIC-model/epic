module interpolation
    use parameters, only : parcel_info, interpl_info
    use parcel_container, only : parcel_container_type, n_parcels
    use ellipse
    use nearest_grid_point
    implicit none

    private :: do_elliptic_to_grid,     &
               do_non_elliptic_to_grid, &
               get_indices_and_weights


    ! interpolation indices
    integer ii(8), jj(8)

    ! interpolation weights
    double precision weight(8)

    private :: ii, jj, weight

    contains

        subroutine to_grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(inout) :: field(:, :, :)

            field = 0.0

            if (parcels_info%is_elliptic) then
                call do_elliptic_to_grid(parcels, attrib, field)
            else
                call do_non_elliptic_to_grid(parcels, attrib, field)
            endif

            ! apply boundary condition

        end subroutine to_grid


        subroutine do_elliptic_to_grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(inout) :: field(:, :, :)
            integer                                    :: ncomp
            double precision                           :: points(2, 2)

            ! number of field components
            ncomp = shape(field)(2)

            do n = 1, n_parcels

                points = get_ellipse_points(parcels%position,   &
                                            parcels%volume,     &
                                            parcels%B11,        &
                                            parcels%B12)

                ! we have 2 points per ellipse
                do point = 1, 2

                    ! get interpolation weights and mesh indices
                    call get_indices(ii, jj, weight)ii

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over field indices
                        do i = 1, len(ii)
                            ! the weight is halved due to 2 points per ellipse
                            field(ii, jj, c) = field(ii, jj, c) + 0.5 * weight(i) * attrib(i)
                        enddo
                    enddo
                enddo
            enddo

        end subroutine do_elliptic_to_grid


        subroutine do_non_elliptic_to_grid(parcels, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(inout) :: field(:, :, :)


        end subroutine do_non_elliptic_to_grid



        subroutine get_indices_and_weights(ii, jj, weight)

            if interpl_info%method == 'nearest-grid-point'

            else if interpl_info% method == 'exact'

            else if interpl_info%method == 'trilinear'

            else
                print *, "Unknown interpolation method '", interpl_info%method, "'."
                stop
            endif

        end subroutine get_indices_and_weights

end module interpolation
