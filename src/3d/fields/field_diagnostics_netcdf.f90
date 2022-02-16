! =============================================================================
!                   Write field diagnostics to NetCDF.
! =============================================================================
module field_diagnostics_netcdf
    use field_diagnostics
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use parameters, only : lower, extent, nx, ny, nz
    use config, only : package_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    implicit none

    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,                   &
                          rms_v_id, abserr_v_id, max_npar_id, min_npar_id, &
                          avg_npar_id, avg_nspar_id

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
            character(:), allocatable :: name

            ncfname =  basename // '_field_stats.nc'

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call read_netcdf_field_stats_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)

            call write_netcdf_attribute(ncid=ncid, name='EPIC_version', val=package_version)
            call write_netcdf_attribute(ncid=ncid, name='file_type', val='field_stats')
            call write_netcdf_attribute(ncid=ncid, name='Conventions', val='CF-1.9')
            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))
            call write_netcdf_timestamp(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='t',                          &
                                         dimsize=NF90_UNLIMITED,            &
                                         dimid=t_dim_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='t',                                                   &
                long_name='time ',                                          &
                std_name='time',                                            &
                unit='seconds since 1970-01-01',                            &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=t_axis_id)

            ncerr = nf90_put_att(ncid, t_axis_id, "axis", 'T')
            call check_netcdf_error("Failed to add axis attribute.")

            ncerr = nf90_put_att(ncid, t_axis_id, "calendar", &
                                 'proleptic_gregorian')
            call check_netcdf_error("Failed to add calendear attribute.")

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='rms_v',                                               &
                long_name='rms volume error',                               &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_v_id)


            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='abserr_v',                                            &
                long_name='max absolute normalised volume error',           &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=abserr_v_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_npar',                                            &
                long_name='max num parcels per cell',                       &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=max_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_npar',                                            &
                long_name='min num parcels per cell',                       &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=min_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_npar',                                            &
                long_name='average num parcels per cell',                   &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_nspar',                                           &
                long_name='average num small parcels per cell',             &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_nspar_id)

            call close_definition(ncid)

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

        end subroutine read_netcdf_field_stats_content

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_stats(t)
            double precision, intent(in)    :: t

            call start_timer(field_stats_io_timer)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, rms_v_id, rms_v, n_writes)
            call write_netcdf_scalar(ncid, abserr_v_id, abserr_v, n_writes)
            call write_netcdf_scalar(ncid, max_npar_id, max_npar, n_writes)
            call write_netcdf_scalar(ncid, min_npar_id, min_npar, n_writes)
            call write_netcdf_scalar(ncid, avg_npar_id, avg_npar, n_writes)
            call write_netcdf_scalar(ncid, avg_nspar_id, avg_nspar, n_writes)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_stats_io_timer)

        end subroutine write_netcdf_field_stats

end module field_diagnostics_netcdf
