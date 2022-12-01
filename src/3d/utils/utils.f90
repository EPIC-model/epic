module utils
    use constants, only : one
    use options, only : field_file          &
                      , field_tol           &
                      , output              &
                      , l_restart           &
                      , restart_file        &
                      , time                &
                      , verbose
    use netcdf_utils, only : set_netcdf_dimensions, set_netcdf_axes
    use field_netcdf
    use field_diagnostics_netcdf
    use field_diagnostics, only : calculate_field_diagnostics
    use parcel_init, only : init_parcels
    use parcel_netcdf
    use parcel_diagnostics_netcdf
    use parcel_diagnostics
    use parcel_container, only : n_parcels, parcel_alloc
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
    use netcdf_reader, only : get_file_type, get_num_steps, get_time, get_netcdf_box
    use parameters, only : lower, extent, update_parameters, read_zeta_boundary_flag, max_num_parcels
    use physics, only : read_physical_quantities, print_physical_quantities, l_peref
    implicit none

    integer :: nfw  = 0    ! number of field writes
    integer :: npw  = 0    ! number of parcel writes
    integer :: nspw = 0    ! number of parcel diagnostics writes
    integer :: nsfw = 0    ! number of field diagnostics writes

    private :: nfw, npw, nspw, nsfw

    contains

        ! Create NetCDF files and set the step number
        subroutine setup_output_files

            if (.not. l_peref) then
                call calculate_peref
            endif

            if (output%write_parcel_stats) then
                call create_netcdf_parcel_stats_file(trim(output%basename), &
                                                     output%overwrite,      &
                                                     l_restart)
            endif

            if (output%write_fields) then
                call create_netcdf_field_file(trim(output%basename), &
                                              output%overwrite,   &
                                              l_restart)
            endif

            if (output%write_field_stats) then
                call create_netcdf_field_stats_file(trim(output%basename),   &
                                                    output%overwrite,        &
                                                    l_restart)
            endif

            if (output%write_parcels) then
                call create_netcdf_parcel_file(trim(output%basename),    &
                                               output%overwrite,         &
                                               l_restart)
            endif

        end subroutine setup_output_files

        ! Write last step to the NetCDF files. For the time step dt, it
        ! writes zero.
        ! @param[in] t is the time
        subroutine write_last_step(t)
            double precision,  intent(in) :: t
            double precision              :: velocity(3, n_parcels)
            double precision              :: strain(5, n_parcels)
            double precision              :: vorticity(3, n_parcels)

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel

            call vorticity_tendency

            call grid2par(velocity, vorticity, strain)

            call calculate_parcel_diagnostics(velocity)
            call calculate_field_diagnostics

            call write_step(t, .true.)
        end subroutine write_last_step

        ! Write step to the NetCDF files.
        ! @param[in] t is the time
        ! @param[in] l_force a logical to force a write (optional)
        subroutine write_step(t, l_force)
            double precision,  intent(in) :: t
            logical, optional, intent(in) :: l_force
            double precision              :: neg = one

            if (present(l_force)) then
                if (l_force) then
                    neg = -one
                endif
            endif

            ! make sure we always write initial setup
            if (output%write_fields .and. &
                (t + epsilon(zero) >= neg * dble(nfw) * output%field_freq)) then
                call write_netcdf_fields(t)
                nfw = nfw + 1
            endif


            if (output%write_parcels .and. &
                (t + epsilon(zero) >= neg * dble(npw) * output%parcel_freq)) then
                call write_netcdf_parcels(t)
                npw = npw + 1
            endif

            if (output%write_parcel_stats .and. &
                (t + epsilon(zero) >= neg * dble(nspw) * output%parcel_stats_freq)) then
                call write_netcdf_parcel_stats(t)
                nspw = nspw + 1
            endif

            if (output%write_field_stats .and. &
                (t + epsilon(zero) >= neg * dble(nsfw) * output%field_stats_freq)) then
                call write_netcdf_field_stats(t)
                nsfw = nsfw + 1
            endif
        end subroutine write_step

        subroutine setup_restart(restart_file, t, file_type)
            character(*),     intent(in)  :: restart_file
            double precision, intent(out) :: t
            character(*),     intent(out) :: file_type
            integer                       :: ncid

            call open_netcdf_file(restart_file, NF90_NOWRITE, ncid)
            call get_file_type(ncid, file_type)
            call get_time(ncid, t)
            call close_netcdf_file(ncid)

            ! set counters (we nee to increment by 1 since
            ! we want to write the next time
            nfw = int(t / output%field_freq) + 1
            npw = int(t / output%parcel_freq) + 1
            nspw = int(t / output%parcel_stats_freq) + 1
            nsfw = int(t / output%field_stats_freq) + 1
        end subroutine setup_restart

        subroutine setup_domain_and_parameters
            character(:), allocatable :: fname
            integer                   :: ncid
            integer                   :: ncells(3)

            ! set axis and dimension names for the NetCDF output
            call set_netcdf_dimensions((/'x', 'y', 'z', 't'/))
            call set_netcdf_axes((/'X', 'Y', 'Z', 'T'/))

            if (l_restart) then
                fname = trim(restart_file)
            else
                fname = trim(field_file)
            endif

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            call get_netcdf_box(ncid, lower, extent, ncells)
            call read_physical_quantities(ncid)
            call read_zeta_boundary_flag(ncid)

            call close_netcdf_file(ncid)

            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)

            ! update global parameters
            call update_parameters

#ifdef ENABLE_VERBOSE
            if (verbose) then
                call print_physical_quantities
            endif
#endif
        end subroutine setup_domain_and_parameters

        subroutine setup_parcels
            character(len=16) :: file_type

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
                time%initial = zero ! make sure user cannot start at arbirtrary time

                call init_parcels(field_file, field_tol)
            endif
        end subroutine setup_parcels

end module utils
