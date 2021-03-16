module interpolation
    use parameters, only : parcel_info, interpl
    use parcel_container, only : parcel_container_type, n_parcels
    use ellipse
    use interpl_methods
    implicit none

    private :: do_elliptic_to_grid,     &
               do_non_elliptic_to_grid, &
               get_indices_and_weights


    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer ij(2, 4)

    ! interpolation weights
    double precision weight(4)

    private :: ij, weight

    contains

        subroutine to_grid(parcels, attrib, field)
            type(parcel_container_type), intent(inout) :: parcels
            double precision,            intent(inout) :: attrib(:, :)
            double precision,            intent(inout) :: field(:, :, :)

            field = 0.0

            if (parcel_info%is_elliptic) then
                call do_elliptic_to_grid(parcels, attrib, field)
            else
                call do_non_elliptic_to_grid(parcels, attrib, field)
            endif

            ! apply boundary condition

        end subroutine to_grid


        subroutine do_elliptic_to_grid(parcels, attrib, field)
            type(parcel_container_type), intent(inout) :: parcels
            double precision,            intent(inout) :: attrib(:, :)
            double precision,            intent(inout) :: field(:, :, :)
            integer                                    :: ncomp, ngp
            double precision                           :: points(2, 2)
            integer                                    :: n, p, c, i
            integer                                    :: the_shape(3)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n),      &
                                            parcels%B11(n),         &
                                            parcels%B12(n))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! get interpolation weights and mesh indices
                    call get_indices_and_weights(points(:, p), ngp)

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            field(ij(1, i), ij(2, i), c) = field(ij(1, i), ij(2, i), c)     &
                                                         + 0.5 * weight(i) * attrib(n, c)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine do_elliptic_to_grid


        subroutine do_non_elliptic_to_grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(inout) :: attrib(:, :)
            double precision,            intent(inout) :: field(:, :, :)


        end subroutine do_non_elliptic_to_grid



        subroutine get_indices_and_weights(pos, ngp)
            double precision, intent(in)    :: pos(2)
            integer,          intent(inout) :: ngp

            if (interpl == 'nearest-grid-point') then
                call nearest_grid_point(pos, ij, weight, ngp)
            else if (interpl == 'trilinear') then
                call trilinear(pos, ij, weight, ngp)
            else
                print *, "Unknown interpolation method '", interpl, "'."
                stop
            endif

        end subroutine get_indices_and_weights

end module interpolation
