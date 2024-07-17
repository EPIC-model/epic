! =============================================================================
!                   Write field diagnostics to NetCDF.
!
! Note: Only the root rank writes field diagnostics.
! =============================================================================
module field_diagnostics_netcdf
    use field_diagnostics
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use constants, only : one
    use parameters, only : lower, extent, nx, ny, nz, write_zeta_boundary_flag
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities, ape_calculation

    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes
    double precision   :: restart_time
    integer            :: field_stats_io_timer

    integer, parameter :: NC_RMS_VOL     = 1    &
                        , NC_ABS_ERR_VOL = 2    &
                        , NC_MAX_NPAR    = 3    &
                        , NC_MIN_NPAR    = 4    &
                        , NC_AVG_NPAR    = 5    &
                        , NC_AVG_NSPAR   = 6    &
                        , NC_KE          = 7    &
                        , NC_EN          = 8    &
                        , NC_APE         = 9    &
                        , NC_MIN_BUOY    = 10   &
                        , NC_MAX_BUOY    = 11

    type(netcdf_info) :: nc_dset(NC_MAX_BUOY)

    public :: create_netcdf_field_stats_file,   &
              write_netcdf_field_stats,         &
              field_stats_io_timer

    contains

        ! Create the NetCDF field diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: n

            if (world%rank .ne. world%root) then
                return
            endif

            call set_netcdf_field_diagnostics_info

            nc_dset(:)%l_enabled = .true.

            nc_dset(NC_APE)%l_enabled = (ape_calculation == 'ape density')

            ncfname =  basename // '_field_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid, l_serial=.true.)
                call get_num_steps(ncid, n_writes)
                if (n_writes > 0) then
                    call get_time(ncid, restart_time)
                    call read_netcdf_field_stats_content
                    call close_netcdf_file(ncid, l_serial=.true.)
                    n_writes = n_writes + 1
                    return
                else
                    call close_netcdf_file(ncid, l_serial=.true.)
                    call delete_netcdf_file(ncfname)
                endif
            endif

            call create_netcdf_file(ncfname, overwrite, ncid, l_serial=.true.)

            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='field_stats',      &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            ! define field diagnostics
            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call define_netcdf_dataset(ncid=ncid,                       &
                                               name=nc_dset(n)%name,            &
                                               long_name=nc_dset(n)%long_name,  &
                                               std_name=nc_dset(n)%std_name,    &
                                               unit=nc_dset(n)%unit,            &
                                               dtype=nc_dset(n)%dtype,          &
                                               dimids=(/t_dim_id/),             &
                                               varid=nc_dset(n)%varid)

                endif
            enddo

            call close_definition(ncid)

            call close_netcdf_file(ncid, l_serial=.true.)

        end subroutine create_netcdf_field_stats_file

        ! Pre-condition: Assumes an open file and nc_dset being initialised.
        subroutine read_netcdf_field_stats_content
            integer :: n

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
                endif
            enddo

        end subroutine read_netcdf_field_stats_content

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_stats(t)
            double precision, intent(in)    :: t

            call start_timer(field_stats_io_timer)

            if (world%rank /= world%root) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            if (t <= restart_time) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid, l_serial=.true.)

            if (n_writes == 1) then
                call write_zeta_boundary_flag(ncid)
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes, l_serial=.true.)

            !
            ! write diagnostics
            !
            call write_diagnostic(NC_RMS_VOL, field_stats(IDX_RMS_V))
            call write_diagnostic(NC_ABS_ERR_VOL, field_stats(IDX_ABSERR_V))
            call write_diagnostic(NC_MAX_NPAR, field_stats(IDX_MAX_NPAR))
            call write_diagnostic(NC_MIN_NPAR, field_stats(IDX_MIN_NPAR))
            call write_diagnostic(NC_AVG_NPAR, field_stats(IDX_AVG_NPAR))
            call write_diagnostic(NC_AVG_NSPAR, field_stats(IDX_AVG_NSPAR))
            call write_diagnostic(NC_KE, field_stats(IDX_KEG))
            call write_diagnostic(NC_APE, field_stats(IDX_APEG))
            call write_diagnostic(NC_EN, field_stats(IDX_ENG))
            call write_diagnostic(NC_MIN_BUOY, field_stats(IDX_MIN_BUOY))
            call write_diagnostic(NC_MAX_BUOY, field_stats(IDX_MAX_BUOY))

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid, l_serial=.true.)

            call stop_timer(field_stats_io_timer)

        end subroutine write_netcdf_field_stats

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine write_diagnostic(id, val)
            integer,          intent(in) :: id
            double precision, intent(in) :: val

            if (nc_dset(id)%l_enabled) then
                if (nc_dset(id)%dtype == NF90_INT) then
                    call write_netcdf_scalar(ncid, nc_dset(id)%varid, int(val), &
                                             n_writes, l_serial=.true.)
                else
                    call write_netcdf_scalar(ncid, nc_dset(id)%varid, val, &
                                             n_writes, l_serial=.true.)
                endif
            endif

        end subroutine write_diagnostic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_diagnostics_info

            nc_dset(NC_RMS_VOL) = netcdf_info(                          &
                name='rms_v',                                           &
                long_name='relative rms volume error',                  &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)


            nc_dset(NC_ABS_ERR_VOL) = netcdf_info(                      &
                name='abserr_v',                                        &
                long_name='max absolute normalised volume error',       &
                std_name='',                                            &
                unit='m^3',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_MAX_NPAR) = netcdf_info(                         &
                name='max_npar',                                        &
                long_name='max num parcels per cell',                   &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_INT)

            nc_dset(NC_MIN_NPAR) = netcdf_info(                         &
                name='min_npar',                                        &
                long_name='min num parcels per cell',                   &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_INT)

            nc_dset(NC_AVG_NPAR) = netcdf_info(                         &
                name='avg_npar',                                        &
                long_name='average num parcels per cell',               &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            nc_dset(NC_AVG_NSPAR) = netcdf_info(                        &
                name='avg_nspar',                                       &
                long_name='average num small parcels per cell',         &
                std_name='',                                            &
                unit='1',                                               &
                dtype=NF90_DOUBLE)

            nc_dset(NC_KE) = netcdf_info(                               &
                name='ke',                                              &
                long_name='domain-averaged kinetic energy',             &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_APE) = netcdf_info(                              &
                name='ape',                                             &
                long_name='domain-averaged available potential energy', &
                std_name='',                                            &
                unit='m^2/s^2',                                         &
                dtype=NF90_DOUBLE)

            nc_dset(NC_EN) = netcdf_info(                               &
                name='en',                                              &
                long_name='domain-averaged enstrophy',                  &
                std_name='',                                            &
                unit='1/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_MIN_BUOY) = netcdf_info(                         &
                name='min_buoyancy',                                    &
                long_name='minimum gridded buoyancy',                   &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

            nc_dset(NC_MAX_BUOY) = netcdf_info(                         &
                name='max_buoyancy',                                    &
                long_name='maximum gridded buoyancy',                   &
                std_name='',                                            &
                unit='m/s^2',                                           &
                dtype=NF90_DOUBLE)

        end subroutine set_netcdf_field_diagnostics_info

end module field_diagnostics_netcdf
