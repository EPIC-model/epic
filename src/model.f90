module model
    use constants, only : max_num_parcels
    use init, only : init_parcels
    use parser, only : read_config_file
    use parcel_container
    use rk4,    only : rk4_step
    use writer, only : h5_create_file
    implicit none

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file

            call alloc_parcel_mem(max_num_parcels)

            call init_parcels

            call h5_create_file

            call h5_write_parcels(0)
        end subroutine


        subroutine run
            use parameters, only : time
            double precision :: t  ! current time
            double precision :: dt ! time step

            dt = time%dt

            do while (t <= time%limit)

                call rk4_step(dt)

                t = t + dt
            end do

        end subroutine run

        subroutine post_run
            call dealloc_parcel_mem
        end subroutine

end module model
