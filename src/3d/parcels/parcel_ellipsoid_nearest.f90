submodule (parcel_ellipsoid) parcel_ellipsoid_nearest_smod
    use mpi_layout, only : box
    use parameters, only : dxi, vmin

    contains

        subroutine parcel_ellipsoid_get_local_cell_index(this, n, ix, iy, iz)
            class(ellipsoid_pc_type), intent(in)  :: this
            integer,                  intent(in)  :: n
            integer,                  intent(out) :: ix, iy, iz

            ix =     int(dxi(1) * (this%position(1, n) - box%halo_lower(1)))
            iy =     int(dxi(2) * (this%position(2, n) - box%halo_lower(2)))
            iz = min(int(dxi(3) * (this%position(3, n) - box%lower(3))), nz-1)

        end subroutine parcel_ellipsoid_get_local_cell_index

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function parcel_ellipsoid_is_small(this, n) result(l_small)
            class(ellipsoid_pc_type), intent(in) :: this
            integer,                  intent(in) :: n
            logical                              :: l_small

            l_small = (this%volume(n) < vmin)

        end function parcel_ellipsoid_is_small

end submodule parcel_ellipsoid_nearest_smod
