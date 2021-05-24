! =============================================================================
!                       EPIC - Elliptic Parcel-in-Cell
! =============================================================================
program epic
    use hdf5
    use constants, only : max_num_parcels, zero
    use diagnostics
    use parser, only : read_config_file, write_h5_params
    use parcel_container
    use parcel_bc
    use parcel_split, only : split_ellipses
    use parcel_merge, only : merge_ellipses
    use parcel_correction, only : init_parcel_correction, apply_laplace, apply_gradient
    use fields
    use tri_inversion, only : init_inversion
    use parcel_interpl
    use rk4
    use model, only : model_init
    use writer, only : open_h5_file,                        &
                       close_h5_file,                       &
                       write_h5_double_scalar_step_attrib,  &
                       write_h5_integer_scalar_step_attrib, &
                       write_h5_integer_scalar_attrib,      &
                       h5err
    implicit none

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

    contains

        subroutine pre_run
            ! parse the config file
            call read_config_file

            call write_h5_params

            call parcel_alloc(max_num_parcels)

            call model_init("TaylorGreen")

            call rk4_alloc(max_num_parcels)

            call init_inversion

            call init_parcel_correction

            call par2grid


        end subroutine


        subroutine run
            use options, only : time, output, verbose, parcel_info
            double precision :: t    = zero ! current time
            double precision :: dt   = zero ! time step
            integer          :: iter = 1    ! simulation iteration
            integer          :: nw   = 0    ! number of writes to h5
            integer          :: cor_iter    ! iterator for parcel correction

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
                    call split_ellipses(parcels, parcel_info%lambda, parcel_info%vmaxfraction)
                endif

                if (mod(iter, parcel_info%correction_freq) == 0) then
                    call vol2grid
                    do cor_iter=1,parcel_info%correction_iters
                        call apply_laplace(volg)
                        call vol2grid
                        if(parcel_info%apply_gradient) then
                            call apply_gradient(volg,parcel_info%gradient_pref)
                            call vol2grid
                        end if
                    end do
                 endif


                t = t + dt
                iter = iter + 1
            end do

            call par2grid

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
            double precision :: H, S11, S12, S21, S22, gmax
            integer          :: i, j

            H = zero
            if (parcel_info%is_elliptic .and. time%is_adaptive) then
                do i = 0, nx-1
                    do j = 0, nz
                        S11 = velgradg(j, i, 1)
                        S12 = velgradg(j, i, 2)
                        S21 = velgradg(j, i, 3)
                        S22 = velgradg(j, i, 4)
                        H = max(H, (S11 - S22) ** 2 + (S12 + S21) ** 2)
                    enddo
                enddo

                gmax = 0.5d0 * dsqrt(H)
                dt = min(time%dt_max, time%alpha / gmax)

            else if (time%is_adaptive) then
                    ! adaptive time stepping according to
                    ! https://doi.org/10.1002/qj.3319
                    max_vorticity = maxval(abs(vortg))
                    dt = min(time%alpha / max_vorticity, time%dt)
            else
                dt = time%dt
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step


    ! Get the file name provided via the command line
    subroutine parse_command_line
        use options, only : filename, verbose
        integer                          :: i
        character(len=32)                :: arg

        filename = ''
        i = 0
        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--config') then
                i = i + 1
                call get_command_argument(i, arg)
                filename = trim(arg)
            else if (arg == '--verbose') then
                verbose = .true.
            endif
            i = i+1
        end do

        if (filename == '') then
            print *, 'No configuration file provided. Run code with "./epic --config [config file]"'
            stop
        endif

        ! This is the main application of EPIC
        if (verbose) then
            print *, 'Running EPIC with "', trim(filename), '"'
        endif
    end subroutine parse_command_line
end program epic
