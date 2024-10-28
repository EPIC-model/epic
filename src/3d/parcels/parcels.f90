module parcels_mod
    use parcel_container, only : pc_type
    use parcel_ellipsoid, only : ellipsoid_pc_type

    type(ellipsoid_pc_type) :: parcels

end module parcels_mod
