submodule (parcel_ellipse) parcel_ellipse_nearest_smod
    use mpi_layout, only : box
    use parameters, only : dxi, amin

    contains

        subroutine parcel_ellipse_get_local_cell_index(this, n, ix, iy, iz)
            class(ellipse_pc_type), intent(in)  :: this
            integer,                intent(in)  :: n
            integer,                intent(out) :: ix, iy, iz

            ix = int(dxi(1) * (this%position(1, n) - box%halo_lower(1)))
            iy = int(dxi(2) * (this%position(2, n) - box%halo_lower(2)))
            iz = 0

        end subroutine parcel_ellipse_get_local_cell_index

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function parcel_ellipse_is_small(this, n) result(l_small)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            logical                            :: l_small

            l_small = (this%area(n) < amin)

        end function parcel_ellipse_is_small

end submodule parcel_ellipse_nearest_smod
