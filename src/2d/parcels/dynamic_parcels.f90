module  dynamic_parcels
use parcel_ellipsoid, only : ellipsoid_parcel_type
implicit none

    class(ellipsoid_parcel_type), allocatable :: parcels
    integer :: n_parcels

end module dynamic_parcels

