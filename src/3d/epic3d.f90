! =============================================================================
!                       EPIC - Elliptical Parcel-in-Cell
! =============================================================================
program epic
    use constants, only : max_num_parcels, zero
    use timer
    use field_diagnostics
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
    use parcel_diagnostics, only : init_parcel_diagnostics,    &
                                   create_h5_parcel_stat_file, &
                                   hdf5_parcel_stat_timer
    use parcel_hdf5
    use fields
    use field_hdf5, only : hdf5_field_timer, create_h5_field_file
    use tri_inversion, only : init_inversion, vor2vel_timer, vtend_timer
    use parcel_interpl, only : grid2par_timer, par2grid_timer
#ifndef NDEBUG
    use parcel_interpl, only : sym_vol2grid_timer
#endif
    use parcel_init, only : init_parcels, init_timer
    use ls_rk4, only : ls_rk4_alloc, ls_rk4_dealloc, ls_rk4_step, rk4_timer
    use h5_utils, only : initialise_hdf5, finalise_hdf5
    use utils, only : write_last_step
#ifdef ENABLE_VERBOSE
    use merge_hdf5, only : create_h5_merger_files
#endif
    implicit none

    integer :: epic_timer

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
            use options, only : field_file, field_tol, output, read_config_file

            call register_timer('epic', epic_timer)
            call register_timer('par2grid', par2grid_timer)
            call register_timer('grid2par', grid2par_timer)
            call register_timer('parcel split', split_timer)
            call register_timer('parcel merge', merge_timer)
            call register_timer('laplace correction', lapl_corr_timer)
            call register_timer('gradient correction', grad_corr_timer)
            call register_timer('parcel init', init_timer)
            call register_timer('parcel hdf5', hdf5_parcel_timer)
            call register_timer('parcel diagnostics hdf5', hdf5_parcel_stat_timer)
            call register_timer('field hdf5', hdf5_field_timer)
            call register_timer('field diagnostics hdf5', hdf5_field_stat_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('parcel push', rk4_timer)
            call register_timer('merge nearest', merge_nearest_timer)
            call register_timer('merge tree resolve', merge_tree_resolve_timer)
#ifndef NDEBUG
            call register_timer('symmetric vol2grid', sym_vol2grid_timer)
#endif

            call start_timer(epic_timer)

            call initialise_hdf5

            ! parse the config file
            call read_config_file

            call parcel_alloc(max_num_parcels)

            call init_parcels(field_file, field_tol)

            call ls_rk4_alloc(max_num_parcels)

            call init_inversion

            call init_parcel_correction

            if (output%h5_write_parcel_stats) then
                call init_parcel_diagnostics
                call create_h5_parcel_stat_file(trim(output%h5_basename), &
                                                output%h5_overwrite)
            endif

            call field_default

            if (output%h5_write_fields) then
                call create_h5_field_file(trim(output%h5_basename), output%h5_overwrite)
            endif

            if (output%h5_write_field_stats) then
                call create_h5_field_stats_file(trim(output%h5_basename), output%h5_overwrite)
            endif

            if (output%h5_write_parcels) then
                call create_h5_parcel_file(trim(output%h5_basename), output%h5_overwrite)
            endif

#ifdef ENABLE_MERGER_DUMP
            call create_h5_merger_files(trim(output%h5_basename), output%h5_overwrite)
#endif
        end subroutine


        subroutine run
            use options, only : time, parcel
#ifdef ENABLE_VERBOSE
            use options, only : verbose
#endif
            double precision :: t    = zero ! current time
            integer          :: cor_iter    ! iterator for parcel correction

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

            ! write final step
            call write_last_step(t)

        end subroutine run

        subroutine post_run
            use options, only : output
            call parcel_dealloc
            call ls_rk4_dealloc
            call finalise_hdf5

            call stop_timer(epic_timer)

            call write_time_to_csv(output%h5_basename)
            call print_timer
        end subroutine


    ! Get the file name provided via the command line
    subroutine parse_command_line
        use options, only : filename
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
                print *, 'Run code with "./epic --config [config file]"'
                stop
#ifdef ENABLE_VERBOSE
            else if (arg == '--verbose') then
                verbose = .true.
#endif
            endif
            i = i+1
        end do

        if (filename == '') then
            print *, 'No configuration file provided. Run code with "./epic --config [config file]"'
            stop
        endif

#ifdef ENABLE_VERBOSE
        ! This is the main application of EPIC
        if (verbose) then
            print *, 'Running EPIC with "', trim(filename), '"'
        endif
#endif
    end subroutine parse_command_line
end program epic
