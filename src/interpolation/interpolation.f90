module interpolation
    use constants, only : max_num_parcels
    use parameters, only : parcel_info, interpl
    use parcel_container, only : parcel_container_type, n_parcels
    use ellipse
    use interpl_methods
    use field_bc
    implicit none

    private :: par2grid_elliptic,       &
               par2grid_non_elliptic,   &
               grid2par_elliptic,       &
               grid2par_non_elliptic,   &
               get_indices_and_weights


    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer ij(2, 4)

    ! interpolation weights
    double precision weight(4)

    private :: ij, weight

    contains

        subroutine par2grid(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(:, :, :)

            field = 0.0

            if (parcel_info%is_elliptic) then
                call par2grid_elliptic(parcels, attrib, field)
            else
                call par2grid_non_elliptic(parcels, attrib, field)
            endif

            ! apply boundary condition
            call apply_field_bc(field)

        end subroutine par2grid

        subroutine par2grid_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
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
                                            parcels%volume(n, 1),   &
                                            parcels%B(n, :))

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
        end subroutine par2grid_elliptic


        subroutine par2grid_non_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(in)    :: parcels
            double precision,            intent(in)    :: attrib(:, :)
            double precision,            intent(inout) :: field(:, :, :)
            integer                                    :: ncomp, ngp
            integer                                    :: n, c, i
            integer                                    :: the_shape(3)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                ! get interpolation weights and mesh indices
                call get_indices_and_weights(parcels%position(n, :), ngp)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        ! the weight is halved due to 2 points per ellipse
                        field(ij(1, i), ij(2, i), c) = field(ij(1, i), ij(2, i), c) &
                                                     + weight(i) * attrib(n, c)
                    enddo
                enddo
            enddo

        end subroutine par2grid_non_elliptic


        subroutine grid2par(parcels, attrib, field)
            type(parcel_container_type), intent(out) :: parcels
            double precision,            intent(out) :: attrib(:, :)
            double precision,            intent(in)  :: field(:, :, :)

            if (parcel_info%is_elliptic) then
                call grid2par_elliptic(parcels, attrib, field)
            else
                call grid2par_non_elliptic(parcels, attrib, field)
            endif

        end subroutine grid2par


        subroutine grid2par_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(out) :: parcels
            double precision,            intent(out) :: attrib(:, :)
            double precision,            intent(in)  :: field(:, :, :)
            integer                                  :: ncomp, ngp
            double precision                         :: points(2, 2)
            integer                                  :: n, p, c, i
            integer                                  :: the_shape(3)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                ! clear old data
                attrib(n, :) = 0.0

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n, 1),   &
                                            parcels%B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! get interpolation weights and mesh indices
                    call get_indices_and_weights(points(:, p), ngp)

                    ! loop over field components
                    do c = 1, ncomp
                        ! loop over grid points which are part of the interpolation
                        do i = 1, ngp
                            ! the weight is halved due to 2 points per ellipse
                            attrib(n, c) = attrib(n, c) &
                                         + 0.5 * weight(i) * field(ij(1, i), ij(2, i), c)
                        enddo
                    enddo
                enddo
            enddo

        end subroutine grid2par_elliptic

        subroutine grid2par_non_elliptic(parcels, attrib, field)
            type(parcel_container_type), intent(out) :: parcels
            double precision,            intent(out) :: attrib(:, :)
            double precision,            intent(in)  :: field(:, :, :)
            integer                                  :: ncomp, ngp
            integer                                  :: n, c, i
            integer                                  :: the_shape(3)

            ! number of field components
            the_shape = shape(field)
            ncomp = the_shape(3)

            do n = 1, n_parcels

                ! clear old data
                attrib(n, :) = 0.0

                ! get interpolation weights and mesh indices
                call get_indices_and_weights(parcels%position(n, :), ngp)

                ! loop over field components
                do c = 1, ncomp
                    ! loop over grid points which are part of the interpolation
                    do i = 1, ngp
                        attrib(n, c) = attrib(n, c) &
                                     + weight(i) * field(ij(1, i), ij(2, i), c)
                    enddo
                enddo
            enddo

        end subroutine grid2par_non_elliptic


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
