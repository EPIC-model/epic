module utils
    use constants, only : one
    use field_io
    use field_diagnostics_io
    use parcel_io
    use parcel_diagnostics_io
    use parcel_diagnostics
    use parcel_container, only : n_parcels
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
    implicit none

    integer :: nfw  = 0    ! number of field writes to h5
    integer :: npw  = 0    ! number of parcel writes to h5
    integer :: nspw = 0    ! number of parcel diagnostics writes to h5
    integer :: nsfw = 0    ! number of field diagnostics writes to h5

    private :: nfw, npw, nspw, nsfw

    contains

        ! Create H5 files and set the step number
        subroutine setup_output_files
            use options, only : output, l_restart

            if (output%write_parcel_stats) then
                call create_parcel_stats_file(trim(output%basename), &
                                              output%overwrite,      &
                                              l_restart)
            endif

            if (output%write_fields) then
                call create_field_file(trim(output%basename), &
                                          output%overwrite,   &
                                          l_restart)
            endif

            if (output%write_field_stats) then
                call create_field_stats_file(trim(output%basename),   &
                                             output%overwrite,        &
                                             l_restart)
            endif

            if (output%write_parcels) then
                call create_parcel_file(trim(output%basename),    &
                                        output%overwrite,         &
                                        l_restart)
            endif

        end subroutine setup_output_files

        ! Write last step to the H5 files. For the time step dt, it
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
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(vortg, tbuoyg, velgradg, vtend)

            call grid2par(velocity, vorticity, strain)

            call calculate_parcel_diagnostics(velocity)

            call write_step(t, .true.)
        end subroutine write_last_step

        ! Write step to the H5 files.
        ! @param[in] t is the time
        ! @param[in] l_force a logical to force a write (optional)
        subroutine write_step(t, l_force)
            use options, only : output
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
                call write_fields(t)
                nfw = nfw + 1
            endif


            if (output%write_parcels .and. &
                (t + epsilon(zero) >= neg * dble(npw) * output%parcel_freq)) then
                call write_parcels(t)
                npw = npw + 1
            endif

            if (output%write_parcel_stats .and. &
                (t + epsilon(zero) >= neg * dble(nspw) * output%parcel_stats_freq)) then
                call write_parcel_stats(t)
                nspw = nspw + 1
            endif

            if (output%write_field_stats .and. &
                (t + epsilon(zero) >= neg * dble(nsfw) * output%field_stats_freq)) then
                call write_field_stats(t)
                nsfw = nsfw + 1
            endif
        end subroutine write_step

        subroutine setup_restart(restart_file, t, file_type)
#ifndef ENABLE_NETCDF
            use hdf5
            use h5_utils, only : initialise_hdf5, finalise_hdf5, open_h5_file, close_h5_file
            use h5_reader, only : get_file_type, get_num_steps, get_time
            integer(hid_t)            :: h5handle
#endif
            character(*),              intent(in)  :: restart_file
            double precision,          intent(out) :: t
            character(:), allocatable, intent(out) :: file_type
            integer                   :: n_steps

#ifndef ENABLE_NETCDF
            call open_h5_file(restart_file, H5F_ACC_RDONLY_F, h5handle)
            call get_file_type(h5handle, file_type)
            call get_num_steps(h5handle, n_steps)
            call get_time(h5handle, n_steps - 1, t)
            call close_h5_file(h5handle)
#endif
        end subroutine setup_restart

end module utils
