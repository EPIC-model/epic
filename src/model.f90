module model
    use constants, only : max_parcel_count
    use init, only : init_parcels
    use parser, only : read_config_file
    use parcels
    use types,  only : time_info_type, parcel_info_type
    use rk4,    only : rk4_step
    implicit none

    type(time_info_type) time_info

    type(parcel_info_type) parcel_info

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file(time_info, &
                                  parcel_info)

            call alloc_parcel_mem(max_parcel_count)

            call init_parcels(parcel_info)
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

        subroutine post_run
            call dealloc_parcel_mem
        end subroutine

end module model
