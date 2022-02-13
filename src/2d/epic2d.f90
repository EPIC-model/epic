! =============================================================================
!                       EPIC2D - Elliptical Parcel-in-Cell
! =============================================================================
program epic2d
    use constants, only : max_num_parcels, zero
    use timer
!     use field_diagnostics, only : hdf5_field_stat_timer
    use parcel_container
    use parcel_bc
    use parcel_split, only : split_ellipses, split_timer
    use parcel_merge, only : merge_ellipses, merge_timer
    use parcel_nearest, only : merge_nearest_timer, merge_tree_resolve_timer
    use parcel_correction, only : init_parcel_correction, &
                                  apply_laplace,          &
                                  apply_gradient,         &
                                  lapl_corr_timer,        &
                                  grad_corr_timer
    use parcel_diagnostics, only : init_parcel_diagnostics, &
                                   parcel_stats_timer
    use parcel_netcdf
    use parcel_diagnostics_netcdf, only : parcel_stats_io_timer
    use fields
    use field_netcdf, only : field_io_timer
    use field_diagnostics, only : field_stats_timer
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use tri_inversion, only : init_inversion, vor2vel_timer, vtend_timer
    use parcel_interpl, only : grid2par_timer, par2grid_timer
#ifndef NDEBUG
    use parcel_interpl, only : sym_vol2grid_timer
#endif
    use parcel_init, only : init_parcels, init_timer
    use ls_rk4, only : ls_rk4_alloc, ls_rk4_dealloc, ls_rk4_step, rk4_timer
    use utils, only : write_last_step, setup_output_files, setup_restart
    implicit none

    integer          :: epic_timer

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
            use options, only : field_file          &
                              , field_tol           &
                              , output              &
                              , read_config_file    &
                              , l_restart           &
                              , restart_file        &
                              , time
            character(16) :: file_type

            call register_timer('epic', epic_timer)
            call register_timer('par2grid', par2grid_timer)
            call register_timer('grid2par', grid2par_timer)
            call register_timer('parcel split', split_timer)
            call register_timer('parcel merge', merge_timer)
            call register_timer('laplace correction', lapl_corr_timer)
            call register_timer('gradient correction', grad_corr_timer)
            call register_timer('parcel initialisation', init_timer)
            call register_timer('parcel diagnostics', parcel_stats_timer)
            call register_timer('parcel I/O', parcel_io_timer)
            call register_timer('parcel diagnostics I/O', parcel_stats_io_timer)
            call register_timer('field I/O', field_io_timer)
            call register_timer('field diagnostics', field_stats_timer)
            call register_timer('field diagnostics I/O', field_stats_io_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('parcel push', rk4_timer)
            call register_timer('merge nearest', merge_nearest_timer)
            call register_timer('merge tree resolve', merge_tree_resolve_timer)
#ifndef NDEBUG
            call register_timer('symmetric vol2grid', sym_vol2grid_timer)
#endif

            call start_timer(epic_timer)

            ! parse the config file
            call read_config_file

            call parcel_alloc(max_num_parcels)

            if (l_restart) then
                call setup_restart(trim(restart_file), time%initial, file_type)

                if (file_type == 'fields') then
                    call init_parcels(restart_file, field_tol)
                else if (file_type == 'parcels') then
                    call read_netcdf_parcels(restart_file)
                else
                    print *, 'Restart file must be of type "fields" or "parcels".'
                    stop
                endif
            else
                time%initial = zero ! make sure user cannot start at arbitrary time

                call init_parcels(field_file, field_tol)
            endif

            call ls_rk4_alloc(max_num_parcels)

            call init_inversion

            call init_parcel_correction

            if (output%write_parcel_stats) then
                call init_parcel_diagnostics
            endif

            call field_default

            call setup_output_files

        end subroutine


        subroutine run
            use options, only : time, parcel
#ifdef ENABLE_VERBOSE
            use options, only : verbose
#endif
            double precision :: t = zero ! current time
            integer          :: cor_iter    ! iterator for parcel correction

            t = time%initial

            do while (t < time%limit)

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print "(a15, f0.4)", "time:          ", t
                endif
#endif
                call ls_rk4_step(t)

                call merge_ellipses(parcels)

                call split_ellipses(parcels, parcel%lambda_max)

                do cor_iter = 1, parcel%correction_iters
                    call apply_laplace
                    call apply_gradient(parcel%gradient_pref, parcel%max_compression)
                enddo

            enddo

            ! write final step (we only write if we really advanced in time)
            if (t > time%initial) then
                call write_last_step(t)
            endif

        end subroutine run

        subroutine post_run
            use options, only : output
            call parcel_dealloc
            call ls_rk4_dealloc

            call stop_timer(epic_timer)

            call write_time_to_csv(output%basename)
            call print_timer
        end subroutine


    ! Get the file name provided via the command line
    subroutine parse_command_line
        use options, only : filename, l_restart, restart_file
#ifdef ENABLE_VERBOSE
        use options, only : verbose
#endif
        integer                          :: i
        character(len=512)               :: arg

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
            else if (arg == '--help') then
                print *, 'Run code with "./epic2d --config [config file]"'
                stop
            else if (arg == '--restart') then
                l_restart = .true.
                i = i + 1
                call get_command_argument(i, arg)
                restart_file = trim(arg)
#ifdef ENABLE_VERBOSE
            else if (arg == '--verbose') then
                verbose = .true.
#endif
            endif
            i = i+1
        end do

        if (filename == '') then
            print *, 'No configuration file provided. Run code with "./epic2d --config [config file]"'
            stop
        endif

        if (l_restart .and. (restart_file == '')) then
            print *, 'No restart file provided. Run code with "./epic2d --config [config file]' // &
                     ' --restart [restart file]"'
            stop
        endif

#ifdef ENABLE_VERBOSE
        ! This is the main application of EPIC
        if (verbose) then
            if (l_restart) then
                print *, 'Restarting EPIC2D with "', trim(filename), '" and "', trim(restart_file), "'"
            else
                print *, 'Running EPIC2D with "', trim(filename), '"'
            endif
        endif
#endif
    end subroutine parse_command_line
end program epic2d
