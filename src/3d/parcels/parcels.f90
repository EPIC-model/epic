module parcels_mod
    use parcel_ellipsoid, only : ellipsoid_pc_type
    use parcel_ellipse, only : ellipse_pc_type

    type(ellipsoid_pc_type) :: parcels      ! interior parcels (3D)
    type(ellipse_pc_type)   :: top_parcels  ! upper surface parcels (2D)
    type(ellipse_pc_type)   :: bot_parcels  ! lower surface parcels (2D)

end module parcels_mod
