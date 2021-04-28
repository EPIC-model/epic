! =============================================================================
!             This module does the setup and runs the simulation.
! =============================================================================
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
    use parcel_diverge, only : init_diverge, apply_diverge
    use fields
    use interpolation
    use rk4
    use writer, only : open_h5_file,                        &
                       close_h5_file,                       &
                       write_h5_double_scalar_step_attrib,  &
                       write_h5_integer_scalar_step_attrib, &
                       write_h5_integer_scalar_attrib,      &
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

            call init_diverge

            ! update volume on the grid
            call par2grid(parcels, parcels%volume, volume_f)

        end subroutine


        subroutine run
            use options, only : time, output, verbose, parcel_info
            double precision :: t    = 0.0 ! current time
            double precision :: dt   = 0.0 ! time step
            integer          :: iter = 1   ! simulation iteration
            integer          :: nw   = 0   ! number of writes to h5

            do while (t <= time%limit)

                dt = get_time_step()

                if (verbose) then
                    print "(a15, f0.4)", "time:          ", t
                    print "(a15, i0)", "iteration:     ", iter
                endif

                ! make sure we always write initial setup
                if (mod(iter - 1, output%h5freq) == 0) then
                    call write_h5_step(nw, t, dt)
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

!                 if (mod(iter, parcel_info%diverge_freq) == 0) then
!                     call apply_diverge(volume_f)
!                 endif

                t = t + dt
                iter = iter + 1
            end do

            ! write final step
            call write_h5_step(nw, t, dt)

        end subroutine run

        subroutine post_run
            call parcel_dealloc
            call rk4_dealloc
        end subroutine


        subroutine write_h5_step(nw, t, dt)
            use options, only : output
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            if (verbose) then
                print "(a30)", "write fields and parcels to h5"
            endif

            call open_h5_file(trim(output%h5fname))

            call write_h5_double_scalar_step_attrib(nw, "t", t)

            call write_h5_double_scalar_step_attrib(nw, "dt", dt)

            call write_h5_integer_scalar_step_attrib(nw, "num parcel", n_parcels)

            call write_h5_parcels(nw)

            call write_h5_diagnostics(nw)

            call write_h5_fields(nw)

            ! update number of iterations to h5 file
            call write_h5_num_steps(nw+1)

            call close_h5_file

            ! increment counter
            nw = nw + 1

        end subroutine write_h5_step

        subroutine write_h5_num_steps(nw)
            use options, only : output
            integer, intent(in) :: nw
            integer(hid_t)      :: group
            logical             :: attr_exists

            group = open_h5_group("/")

            ! in principle not necessary but we still check
            call h5aexists_f(group, "nsteps", attr_exists, h5err)

            if (attr_exists) then
                call h5adelete_f(group, "nsteps", h5err)
            endif

            call write_h5_integer_scalar_attrib(group, "nsteps", nw)

            call h5gclose_f(group, h5err)

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
