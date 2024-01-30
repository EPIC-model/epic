module parcels_mod
    use parcel_container, only : pc_type
    use parcel_ellipsoid, only : ellipsoid_pc_type
    use parcel_ellipse, only : ellipse_pc_type
    use parcel_merging, only : parcel_ellipsoid_merge
    use parcel_ellipse_merge, only : surface_parcel_merge

    type(ellipsoid_pc_type) :: parcels      ! interior parcels (3D)
    type(ellipse_pc_type)   :: top_parcels  ! upper surface parcels (2D)
    type(ellipse_pc_type)   :: bot_parcels  ! lower surface parcels (2D)

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_merge

            call parcel_ellipsoid_merge

            call surface_parcel_merge(bot_parcels)
            call surface_parcel_merge(top_parcels)

        end subroutine parcel_merge

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


end module parcels_mod
