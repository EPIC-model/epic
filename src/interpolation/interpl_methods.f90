module interpl_methods
    use parameters, only : mesh
    use fields, only : get_mesh_spacing, get_lower_index
    implicit none

    contains

        !
        ! nearest grid point (NGP) interpolation
        !
        subroutine nearest_grid_point(pos, idx, weight)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: idx(2)
            double precision, intent(out) :: weight
            double precision              :: xy(2)
            integer                       :: i

            idx = get_lower_index(pos)

            dx = get_mesh_spacing()

            xy = mesh%origin + idx * mesh

            do i = 1, 2
                if (pos(i) - xy(i) > 0.5 * dx(i)) then
                    idx(i) = idx(i) + 1
                endif
            enddo

            weight = 1.0

        end subroutine nearest_grid_point

        !
        ! exact "interpolation"
        !
        subroutine exact(pos)
            double precision, intent(inout)  :: pos(2)
            integer,          intent(out)    :: weight

        end subroutine exact

end interpl_methods
