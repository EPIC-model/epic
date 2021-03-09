module model
    use parser, only : read_config_file
    use types, only : time_info_type, parcel_info_type
    implicit none

    type(time_info_type) time_info

    type(parcel_info_type) parcel_info

    contains
        subroutine setup
            ! parse the config file
            call read_config_file(time_info, parcel_info)

        end subroutine


        subroutine run

        end subroutine run


end module model
