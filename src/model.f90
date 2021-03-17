module model
    use constants, only : max_num_parcels
    use init, only : init_parcels, init_velocity_field
    use parser, only : read_config_file, write_h5_params
    use parcel_container
    use parcel_bc
    use fields
    use interpolation
    use rk4
    use writer, only : open_h5_file,             &
                       close_h5_file,            &
                       write_h5_scalar_attrib
    implicit none

    private :: write_h5_step

    contains
        subroutine pre_run
            ! parse the config file
            call read_config_file

            call write_h5_params

            call parcel_alloc(max_num_parcels)

            call rk4_alloc(max_num_parcels)

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
                    call write_h5_step(iter, t, dt)
                endif

                call rk4_step(dt)

                t = t + dt
                iter = iter + 1
            end do

        end subroutine run

        subroutine post_run
            call parcel_dealloc
            call rk4_dealloc
        end subroutine


        subroutine write_h5_step(iter, t, dt)
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

        end subroutine write_h5_step

end module model
