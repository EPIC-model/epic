module  dynamic_parcels
use parcel_types, only : idealised_parcel_type
implicit none

    type(idealised_parcel_type) :: parcels
    integer :: n_parcels

end module dynamic_parcels

