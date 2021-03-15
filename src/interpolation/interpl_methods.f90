module interpl_methods
    use parameters, only : mesh
    use fields, only : get_mesh_spacing, get_lower_index
    implicit none

    contains

        !
        ! nearest grid point (NGP) interpolation
        !
        subroutine nearest_grid_point(pos, ij, weight, ngp)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ij(8, 2)
            double precision, intent(out) :: weight(8)
            integer,          intent(out) :: ngp
            double precision              :: xy(2), dx(2)
            integer                       :: i

            ij(1, :) = get_lower_index(pos)

            dx = get_mesh_spacing()

            xy = mesh%origin + ij(1, :) * dx

            do i = 1, 2
                if (pos(i) - xy(i) > 0.5 * dx(i)) then
                    ij(1, i) = ij(1, i) + 1
                endif
            enddo

            weight(1) = 1.0
            ngp = 1

        end subroutine nearest_grid_point

        !
        ! exact "interpolation"
        !
!         subroutine exact(pos)
!             double precision, intent(inout)  :: pos(2)
!             integer,          intent(out)    :: weight
!
!         end subroutine exact

end module interpl_methods
