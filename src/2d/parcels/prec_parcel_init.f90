! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module prec_parcel_init
    use precipitation_parcels, only : prec_parcels
    use parcel_types, only : prec_parcel_type
    use parameters, only : max_num_prec_parcels

    use omp_lib
    implicit none

    contains

        subroutine initiate_prec_parcel_type

            if (allocated(prec_parcels)) deallocate(prec_parcels)
            allocate(prec_parcel_type :: prec_parcels)

            prec_parcels%dim_string='xz'
            call prec_parcels%set_dimensions()
            call prec_parcels%alloc(max_num_prec_parcels)

        end subroutine initiate_prec_parcel_type

end module prec_parcel_init
