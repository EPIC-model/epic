module model
    use hdf5
    use constants, only : max_num_parcels
    use diagnostics
    use init, only : init_parcels, init_fields
    use parser, only : read_config_file, write_h5_params
    use parcel_container
    use parcel_bc
    use parcel_split, only : split_ellipses
    use parcel_merge, only : merge_ellipses
    use fields
    use interpolation
    use rk4
    use writer, only : open_h5_file,                   &
                       close_h5_file,                  &
                       write_h5_scalar_step_attrib,    &
                       write_h5_integer_scalar_attrib, &
                       h5err
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

            call init_fields

            ! update volume on the grid
            call par2grid(parcels, parcels%volume, volume_f)

        end subroutine


        subroutine run
            use options, only : time, output, verbose, parcel_info
            double precision :: t    = 0.0 ! current time
            double precision :: dt   = 0.0 ! time step
            integer          :: iter = 0   ! simulation iteration
            integer          :: nw   = 0   ! number of writes to h5

            do while (t <= time%limit)

                dt = get_time_step()

                if (verbose) then
                    print "(a15, f0.4)", "time:          ", t
                    print "(a15, i0)", "iteration:     ", iter
                endif

                if (mod(iter, output%h5freq) == 0) then
                    if(verbose) then
                        print "(a30)", "write fields and parcels to h5"
                    endif
                    call write_h5_step(nw, t, dt)
                    nw = nw + 1
                endif

                call rk4_step(dt)

                if (parcel_info%is_elliptic .and.           &
                    mod(iter, parcel_info%merge_freq) == 0) then
                    call merge_ellipses(parcels)
                endif

                if (parcel_info%is_elliptic .and.           &
                    mod(iter, parcel_info%split_freq) == 0) then
                    call split_ellipses(parcels, parcel_info%lambda)
                endif


                ! update volume on the grid
                call par2grid(parcels, parcels%volume, volume_f)

                t = t + dt
                iter = iter + 1
            end do

            if (mod(iter, output%h5freq) == 0) then
                if(verbose) then
                    print "(a30)", "write fields and parcels to h5"
                endif
                call write_h5_step(nw, t, dt)
                nw = nw + 1
            endif

            ! write number of iterations to h5 file
            call write_h5_num_steps(nw)

        end subroutine run

        subroutine post_run
            call parcel_dealloc
            call rk4_dealloc
        end subroutine


        subroutine write_h5_step(iter, t, dt)
            use options, only : output
            integer,          intent(in) :: iter
            double precision, intent(in) :: t
            double precision, intent(in) :: dt

            call open_h5_file(trim(output%h5fname))

            call write_h5_scalar_step_attrib(iter, "t", t)

            call write_h5_scalar_step_attrib(iter, "dt", dt)

            call write_h5_parcels(iter)

            call write_h5_diagnostics(iter)

            call write_h5_fields(iter)

            call close_h5_file

        end subroutine write_h5_step

        subroutine write_h5_num_steps(iter)
            use options, only : output
            integer, intent(in) :: iter
            integer(hid_t)      :: group

            call open_h5_file(trim(output%h5fname))

            group = open_h5_group("/")

            call write_h5_integer_scalar_attrib(group, "nsteps", iter)

            call h5gclose_f(group, h5err)

            call close_h5_file

        end subroutine write_h5_num_steps


        function get_time_step() result(dt)
            use options, only : time
            double precision :: dt
            double precision :: max_vorticity

            if (time%is_adaptive) then
                ! adaptive time stepping according to
                ! https://doi.org/10.1002/qj.3319
                max_vorticity = maxval(abs(vorticity_f))
                dt = min(0.5 / max_vorticity, time%dt)
            else
                dt = time%dt
            endif

            if (dt <= 0.0) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step

end module model
