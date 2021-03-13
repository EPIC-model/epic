module model
    use constants, only : max_num_parcels
    use init, only : init_parcels, init_velocity_field
    use parser, only : read_config_file
    use parcel_container
    use fields
    use rk4,    only : rk4_step
    use writer, only : open_h5_file,             &
                       close_h5_file,            &
                       write_h5_scalar_attrib
    implicit none

    private :: write_step_to_h5

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file

            call alloc_parcel_mem(max_num_parcels)

            call init_parcels

            call init_velocity_field

        end subroutine


        subroutine run
            use parameters, only : time, output
            double precision :: t    = 0.0 ! current time
            double precision :: dt   = 0.0 ! time step
            integer          :: iter = 0

            dt = time%dt

            do while (t <= time%limit)

                if (mod(iter, output%h5freq) == 0) then
                    call write_step_to_h5(iter, t, dt)
                endif

                call rk4_step(dt)

                t = t + dt
                iter = iter + 1
            end do

        end subroutine run

        subroutine post_run
            call dealloc_parcel_mem
        end subroutine


        subroutine write_step_to_h5(iter, t, dt)
            use parameters, only : output
            integer,          intent(in) :: iter
            double precision, intent(in) :: t
            double precision, intent(in) :: dt

            call open_h5_file(trim(output%h5fname))

            call write_h5_scalar_attrib(iter, "t", t)

            call write_h5_scalar_attrib(iter, "dt", dt)

            call write_h5_parcels(iter)

            call write_h5_fields(iter)

            call close_h5_file

        end subroutine write_step_to_h5

end module model
