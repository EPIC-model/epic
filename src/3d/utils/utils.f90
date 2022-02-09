module utils
    use constants, only : one
    use field_hdf5
    use parcel_hdf5
    use parcel_diagnostics
    use field_diagnostics
    use parcel_container, only : n_parcels
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
    use field_diagnostics, only : write_h5_field_stats_step
#ifdef ENABLE_NETCDF
    use field_netcdf
#endif
    implicit none

    integer :: nfw  = 0    ! number of field writes to h5
    integer :: nnfw = 1
    integer :: npw  = 0    ! number of parcel writes to h5
    integer :: nspw = 0    ! number of parcel diagnostics writes to h5
    integer :: nsfw = 0    ! number of field diagnostics writes to h5

    private :: nfw, npw, nspw, nsfw, nnfw

    contains

        ! Create H5 files and set the step number
        subroutine setup_output_files
            use options, only : output, l_restart

            if (output%h5_write_parcel_stats) then
                call create_h5_parcel_stat_file(trim(output%h5_basename), &
                                                output%h5_overwrite,      &
                                                l_restart, nspw)
            endif

            if (output%h5_write_fields) then
                call create_h5_field_file(trim(output%h5_basename), &
                                          output%h5_overwrite,      &
                                          l_restart, nfw)
#ifdef ENABLE_NETCDF
                call create_netcdf_field_file(trim(output%h5_basename), output%h5_overwrite)
#endif
            endif

            if (output%h5_write_field_stats) then
                call create_h5_field_stats_file(trim(output%h5_basename),   &
                                                output%h5_overwrite,        &
                                                l_restart, nsfw)
            endif

            if (output%h5_write_parcels) then
                call create_h5_parcel_file(trim(output%h5_basename),    &
                                           output%h5_overwrite,         &
                                           l_restart, npw)
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

            call calc_parcel_diagnostics(velocity)

            call write_step(t, zero, .true.)
        end subroutine write_last_step

        ! Write step to the H5 files.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        ! @param[in] l_force a logical to force a write (optional)
        subroutine write_step(t, dt, l_force)
            use options, only : output
            double precision,  intent(in) :: t
            double precision,  intent(in) :: dt
            logical, optional, intent(in) :: l_force
            double precision              :: neg = one

            if (present(l_force)) then
                if (l_force) then
                    neg = -one
                endif
            endif

            ! make sure we always write initial setup
            if (output%h5_write_fields .and. &
                (t + epsilon(zero) >= neg * dble(nfw) * output%h5_field_freq)) then
                call write_h5_field_step(nfw, t, dt)
                call write_netcdf_field_step(nnfw, t, dt)
            endif


            if (output%h5_write_parcels .and. &
                (t + epsilon(zero) >= neg * dble(npw) * output%h5_parcel_freq)) then
                call write_h5_parcel_step(npw, t, dt)
            endif

            if (output%h5_write_parcel_stats .and. &
                (t + epsilon(zero) >= neg * dble(nspw) * output%h5_parcel_stats_freq)) then
                call write_h5_parcel_stats_step(nspw, t, dt)
            endif

            if (output%h5_write_field_stats .and. &
                (t + epsilon(zero) >= neg * dble(nsfw) * output%h5_field_stats_freq)) then

                call write_h5_field_stats_step(nsfw, t, dt)
            endif
        end subroutine write_step

end module utils
