module model
    use constants, only : max_num_parcels
    use init, only : init_parcels
    use parser, only : read_config_file
    use parcel_container
    use rk4,    only : rk4_step
    use writer, only : open_h5_file,             &
                       close_h5_file,            &
                       write_h5_scalar_attrib
    implicit none

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file

            call alloc_parcel_mem(max_num_parcels)

            call init_parcels

        end subroutine


        subroutine run
            use parameters, only : time
            double precision :: t    = 0.0 ! current time
            double precision :: dt   = 0.0 ! time step
            integer          :: iter = 0

            dt = time%dt

            do while (t <= time%limit)

                call open_h5_file("test.h5")

                call write_h5_scalar_attrib(iter, "t", t)

                call write_h5_scalar_attrib(iter, "dt", dt)

                call write_h5_parcels(iter)

                call close_h5_file

                call rk4_step(dt)

                t = t + dt
                iter = iter + 1
            end do

        end subroutine run

        subroutine post_run
            call dealloc_parcel_mem
        end subroutine

end module model
