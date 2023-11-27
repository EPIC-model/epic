module utils
    use constants, only : zero, one
    use options, only : field_file          &
                      , flux_file           &
                      , output              &
                      , l_restart           &
                      , restart_file        &
                      , time                &
                      , verbose
    use netcdf_utils, only : set_netcdf_dimensions, set_netcdf_axes
    use field_netcdf
    use fields, only : vortg, field_default, tbuoyg
    use field_diagnostics_netcdf
    use field_diagnostics, only : calculate_field_diagnostics
    use parcel_init, only : parcel_default, init_parcels_from_grids
    use parcel_netcdf
    use parcel_diagnostics_netcdf
    use parcel_diagnostics
    use parcel_container, only : n_parcels
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
    use netcdf_reader, only : get_file_type, get_num_steps, get_time, get_netcdf_box
    use parameters, only : lower, extent, update_parameters, read_zeta_boundary_flag &
                         , set_zeta_boundary_flag
    use bndry_fluxes, only : read_bndry_fluxes
    use physics, only : read_physical_quantities        &
#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
                      , calculate_basic_reference_state &
#endif
                      , print_physical_quantities
    use mpi_layout, only : mpi_layout_init
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    double precision :: t_fw  = zero    ! next intended time of field write
    double precision :: t_pw  = zero    ! next intended time of parcel write
    double precision :: t_pdw = zero    ! next intended time of parcel diagnostics write
    double precision :: t_fdw = zero    ! next intended time of field diagnostics write

    private :: t_fw, t_pw, t_pdw, t_fdw, get_next_write_time

    contains

        ! Get the next intended time of write
        function get_next_write_time(freq, t) result(tnext)
            double precision, intent(in) :: freq, t
            double precision             :: tnext

            tnext = freq * dble(int(t / freq) + 1)

        end function get_next_write_time

        ! Create NetCDF files and set the step number
        subroutine setup_output_files

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

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk substep
            call vor2vel

            call vorticity_tendency

            call grid2par

            call calculate_parcel_diagnostics
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
            if (output%write_fields .and. (t + epsilon(zero) >= neg * t_fw)) then
                call write_netcdf_fields(t)
                t_fw = get_next_write_time(output%field_freq, t)
            endif

            if (output%write_parcels .and. (t + epsilon(zero) >= neg * t_pw)) then
                call write_netcdf_parcels(t)
                t_pw = get_next_write_time(output%parcel_freq, t)
            endif

            if (output%write_parcel_stats .and. (t + epsilon(zero) >= neg * t_pdw)) then
                call write_netcdf_parcel_stats(t)
                t_pdw = get_next_write_time(output%parcel_stats_freq, t)
            endif

            if (output%write_field_stats .and. (t + epsilon(zero) >= neg * t_fdw)) then
                call write_netcdf_field_stats(t)
                t_fdw = get_next_write_time(output%field_stats_freq, t)
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

            ! set the next intended time to write to files:
            t_fw = get_next_write_time(output%field_freq, t)
            t_pw = get_next_write_time(output%parcel_freq, t)
            t_pdw = get_next_write_time(output%parcel_stats_freq, t)
            t_fdw = get_next_write_time(output%field_stats_freq, t)
        end subroutine setup_restart

        subroutine setup_domain_and_parameters
            character(512) :: fname
            integer        :: ncid
            integer        :: ncells(3)

            ! set axis and dimension names for the NetCDF output
            call set_netcdf_dimensions((/'x', 'y', 'z', 't'/))
            call set_netcdf_axes((/'X', 'Y', 'Z', 'T'/))

            if (l_restart) then
                fname = restart_file
            else
                fname = field_file
            endif

            call open_netcdf_file(trim(fname), NF90_NOWRITE, ncid)

            call get_netcdf_box(ncid, lower, extent, ncells)
            call read_physical_quantities(ncid)
            call read_zeta_boundary_flag(ncid)

            call close_netcdf_file(ncid)

            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)

            call mpi_layout_init(lower, extent, nx, ny, nz)

            ! update global parameters
            call update_parameters

#ifdef ENABLE_VERBOSE
            if (verbose) then
                call print_physical_quantities
            endif
#endif
        end subroutine setup_domain_and_parameters

        ! Reads always the last time step of a field file.
        subroutine setup_fields_and_parcels
            character(len=16) :: file_type
            integer           :: ncid

            call field_default

            call parcel_default

            if (l_restart) then
                call setup_restart(trim(restart_file), time%initial, file_type)

                if (file_type == 'fields') then
                    call read_netcdf_fields(trim(restart_file), -1)
                    call init_parcels_from_grids
                else if (file_type == 'parcels') then
                    call read_netcdf_parcels(restart_file)
                else
                    call mpi_exit_on_error('Restart file must be of type "fields" or "parcels".')
                endif
            else
                time%initial = zero ! make sure user cannot start at arbirtrary time

                call open_netcdf_file(field_file, NF90_NOWRITE, ncid)
                call get_file_type(ncid, file_type)
                call close_netcdf_file(ncid)

                if (file_type == 'fields') then
                    call read_netcdf_fields(field_file, -1)
                    call init_parcels_from_grids

                   ! we must check if zeta must be kept zero
                   ! on a vertical boundary
                   call set_zeta_boundary_flag(vortg(:, :, :, I_Z))
               else if (file_type == 'parcels') then
                   call read_netcdf_parcels(field_file)
               else
                   call mpi_exit_on_error('Input file must be of type "fields" or "parcels".')
               endif
            endif

            call read_bndry_fluxes(trim(flux_file))
            
#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
            ! Calculate the squared buoyancy frequency if not provided.
            call calculate_basic_reference_state(nx, ny, nz, extent(3), tbuoyg)
#endif

        end subroutine setup_fields_and_parcels

end module utils
