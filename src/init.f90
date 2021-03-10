module init
    use parcels
    use types,  only : parcel_info_type
    implicit none

    private :: init_random_positions, &
               init_regular_positions

    contains

        subroutine init_parcels(parcel_info)
            type(parcel_info_type), intent(in) :: parcel_info

            if (parcel_info%is_random) then
                call init_random_positions(parcel_info)
            else
                call init_regular_positions(parcel_info)
            endif
        end subroutine init_parcels


        subroutine init_random_positions(parcel_info)
            type(parcel_info_type), intent(in) :: parcel_info

        end subroutine init_random_positions


        subroutine init_regular_positions(parcel_info)
            type(parcel_info_type), intent(in) :: parcel_info

        end subroutine init_regular_positions

end module init
