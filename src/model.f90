module model
    use constants, only : max_num_parcels
    use init, only : init_parcels
    use parser, only : read_config_file
    use parcel_container
    use rk4,    only : rk4_step
    implicit none

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file

            call alloc_parcel_mem(max_num_parcels)

            call init_parcels
        end subroutine


        subroutine run
            use parameters, only : tmax, dt
            double precision :: t  ! current time

            do while (t <= tmax)

                call rk4_step(dt)

                t = t + dt
            end do

        end subroutine run

        subroutine post_run
            call dealloc_parcel_mem
        end subroutine

end module model
