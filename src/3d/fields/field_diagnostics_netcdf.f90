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
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,                   &
                          rms_v_id, abserr_v_id, max_npar_id, min_npar_id, &
                          avg_npar_id, avg_nspar_id, keg_id, eng_id,       &
                          max_buoy_id, min_buoy_id

    double precision   :: restart_time

    integer :: field_stats_io_timer

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

            if (mpi_rank .ne. mpi_master) then
                return
            endif

            ncfname =  basename // '_field_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid, l_single=.true.)
                call get_num_steps(ncid, n_writes)
                call get_time(ncid, restart_time)
                call read_netcdf_field_stats_content
                call close_netcdf_file(ncid, l_single=.true.)
                n_writes = n_writes + 1
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid, l_single=.true.)

            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='field_stats',      &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='rms_v',                                               &
                long_name='relative rms volume error',                      &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_v_id)


            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='abserr_v',                                            &
                long_name='max absolute normalised volume error',           &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=abserr_v_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_npar',                                            &
                long_name='max num parcels per cell',                       &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=max_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_npar',                                            &
                long_name='min num parcels per cell',                       &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=min_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_npar',                                            &
                long_name='average num parcels per cell',                   &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_nspar',                                           &
                long_name='average num small parcels per cell',             &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_nspar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='kinetic energy',                                 &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=keg_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='en',                                                  &
                long_name='enstrophy',                                      &
                std_name='',                                                &
                unit='m^3/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=eng_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_buoyancy',                                        &
                long_name='minimum gridded buoyancy',                       &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=min_buoy_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_buoyancy',                                        &
                long_name='maximum gridded buoyancy',                       &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=max_buoy_id)

            call close_definition(ncid)

            call close_netcdf_file(ncid, l_single=.true.)

        end subroutine create_netcdf_field_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'rms_v', rms_v_id)

            call get_var_id(ncid, 'abserr_v', abserr_v_id)

            call get_var_id(ncid, 'max_npar', max_npar_id)

            call get_var_id(ncid, 'min_npar', min_npar_id)

            call get_var_id(ncid, 'avg_npar', avg_npar_id)

            call get_var_id(ncid, 'avg_nspar', avg_nspar_id)

            call get_var_id(ncid, 'ke', keg_id)

            call get_var_id(ncid, 'en', eng_id)

            call get_var_id(ncid, 'min_buoyancy', min_buoy_id)

            call get_var_id(ncid, 'max_buoyancy', max_buoy_id)

        end subroutine read_netcdf_field_stats_content

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_stats(t)
            double precision, intent(in)    :: t

            call start_timer(field_stats_io_timer)

            if (mpi_rank .ne. mpi_master) then
                return
            endif

            if (t <= restart_time) then
                call stop_timer(field_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid, l_single=.true.)

            if (n_writes == 1) then
                call write_zeta_boundary_flag(ncid)
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, rms_v_id, field_stats(IDX_RMS_V), n_writes)
            call write_netcdf_scalar(ncid, abserr_v_id, field_stats(IDX_ABSERR_V), n_writes)
            call write_netcdf_scalar(ncid, max_npar_id, field_stats(IDX_MAX_NPAR), n_writes)
            call write_netcdf_scalar(ncid, min_npar_id, field_stats(IDX_MIN_NPAR), n_writes)
            call write_netcdf_scalar(ncid, avg_npar_id, field_stats(IDX_AVG_NPAR), n_writes)
            call write_netcdf_scalar(ncid, avg_nspar_id, field_stats(IDX_AVG_NSPAR), n_writes)
            call write_netcdf_scalar(ncid, keg_id, field_stats(IDX_KEG), n_writes)
            call write_netcdf_scalar(ncid, eng_id, field_stats(IDX_ENG), n_writes)
            call write_netcdf_scalar(ncid, min_buoy_id, field_stats(IDX_MIN_BUOY), n_writes)
            call write_netcdf_scalar(ncid, max_buoy_id, field_stats(IDX_MAX_BUOY), n_writes)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid, l_single=.true.)

            call stop_timer(field_stats_io_timer)

        end subroutine write_netcdf_field_stats

end module field_diagnostics_netcdf
