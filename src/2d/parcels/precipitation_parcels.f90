module  precipitation_parcels
use parcel_types, only : prec_parcel_type
implicit none

    class(prec_parcel_type), allocatable :: prec_parcels
    integer :: n_prec_parcels

end module precipitation_parcels

