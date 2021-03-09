module model
    use parser, only : read_config_file
    use types,  only : time_info_type, parcel_info_type
    use rk4,    only : rk4_step
    implicit none

    type(time_info_type) time_info

    type(parcel_info_type) parcel_info

    contains
        subroutine setup
            ! parse the config file
            call read_config_file(time_info, &
                                  parcel_info)

        end subroutine


        subroutine run
            double precision :: t  ! current time
            double precision :: dt ! time step

            dt = time_info%dt

            do while (t <= time_info%limit)

                call rk4_step(dt)

                t = t + dt
            end do

        end subroutine run


end module model
