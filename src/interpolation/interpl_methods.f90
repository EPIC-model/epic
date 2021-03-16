module interpl_methods
    use parameters, only : mesh
    use fields, only : get_mesh_spacing, &
                       get_lower_index,  &
                       get_position
    implicit none

    contains

        !
        ! nearest grid point (NGP) interpolation
        !
        subroutine nearest_grid_point(pos, ij, weight, ngp)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ij(4, 2)
            double precision, intent(out) :: weight(4)
            integer,          intent(out) :: ngp
            double precision              :: xy(2), dx(2)
            integer                       :: i

            ij(1, :) = get_lower_index(pos)

            dx = get_mesh_spacing()

            xy = get_position(ij(1, :))

            do i = 1, 2
                if (pos(i) - xy(i) > 0.5 * dx(i)) then
                    ij(1, i) = ij(1, i) + 1
                endif
            enddo

            weight(1) = 1.0
            ngp = 1

        end subroutine nearest_grid_point

        !
        ! tri-linear interpolation
        !
        subroutine trilinear(pos, ij, weight, ngp)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ij(4, 2)
            double precision, intent(out) :: weight(4)
            integer,          intent(out) :: ngp
            double precision              :: xy(2), dx(2)

            dx = get_mesh_spacing()

            ! (i, j)
            ij(1, :) = get_lower_index(pos)
            xy = get_position(ij(1, :))
            weight(1) = product(1.0 - abs(pos - xy) / dx)

            ! (i+1, j)
            ij(2, :) = ij(1, :)
            ij(2, 1) = ij(2, 1) + 1
            xy = get_position(ij(2, :))
            weight(2) = product(1.0 - abs(pos - xy) / dx)

            ! (i, j+1)
            ij(3, :) = ij(1, :)
            ij(3, 2) = ij(3, 2) + 1
            xy = get_position(ij(3, :))
            weight(3) = product(1.0 - abs(pos - xy) / dx)

            ! (i+1, j+1)
            ij(4, :) = ij(2, :)
            ij(4, 2) = ij(4, 2) + 1
            xy = get_position(ij(4, :))
            weight(4) = product(1.0 - abs(pos - xy) / dx)

            ngp = 4

        end subroutine trilinear

end module interpl_methods
