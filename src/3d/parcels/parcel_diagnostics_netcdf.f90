! =============================================================================
!                      Write parcel diagnostics to NetCDF
! =============================================================================
module parcel_diagnostics_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics
    use parameters, only : lower, extent, nx, ny, nz
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none


    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          pe_id, ke_id, te_id, npar_id, nspar_id,   &
                          rms_x_vor_id, rms_y_vor_id, rms_z_vor_id, &
                          avg_lam_id, std_lam_id,                   &
                          avg_vol_id, std_vol_id, sum_vol_id,       &
                          psi_id, n_par_split_id, n_par_merge_id

    double precision   :: restart_time

    integer :: parcel_stats_io_timer

    public :: create_netcdf_parcel_stats_file,  &
              write_netcdf_parcel_stats,        &
              parcel_stats_io_timer


    contains

        ! Create the parcel diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            ncfname =  basename // '_parcel_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call get_time(ncid, restart_time)
                call read_netcdf_parcel_stats_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   epic_version=package_version, &
                                   file_type='parcel_stats',     &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='pe',                                                  &
                long_name='potential energy',                               &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=pe_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='kinetic energy',                                 &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=ke_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='te',                                                  &
                long_name='total energy',                                   &
                std_name='',                                                &
                unit='m^5/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=te_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='psi',                                                 &
                long_name='enstrophy',                                      &
                std_name='',                                                &
                unit='m^3/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=psi_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='n_parcels',                                           &
                long_name='number of parcels',                              &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='n_small_parcel',                                      &
                long_name='number of small parcels',                        &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=nspar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_lam',                                             &
                long_name='average aspect ratio',                           &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_lam_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='std_lam',                                             &
                long_name='standard deviation aspect ratio',                &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=std_lam_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg_vol',                                             &
                long_name='average volume',                                 &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_vol_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='std_vol',                                             &
                long_name='standard deviation volume',                      &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=std_vol_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='sum_vol',                                             &
                long_name='total volume',                                   &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=sum_vol_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='x_rms_vorticity',                                     &
                long_name='root mean square of x vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_x_vor_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='y_rms_vorticity',                                     &
                long_name='root mean square of y vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_y_vor_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z_rms_vorticity',                                     &
                long_name='root mean square of z vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_z_vor_id)

            call define_netcdf_dataset(                                      &
                 ncid=ncid,                                                  &
                 name='n_parcel_splits',                                     &
                 long_name='number of parcel splits since last time',        &
                 std_name='',                                                &
                 unit='1',                                                   &
                 dtype=NF90_DOUBLE,                                          &
                 dimids=(/t_dim_id/),                                        &
                 varid=n_par_split_id)
            
            call define_netcdf_dataset(                                      &
                 ncid=ncid,                                                  &
                 name='n_parcel_merges',                                     &
                 long_name='number of parcel merges since last time',        &
                 std_name='',                                                &
                 unit='1',                                                   &
                 dtype=NF90_DOUBLE,                                          &
                 dimids=(/t_dim_id/),                                        &
                 varid=n_par_merge_id)
            
            call close_definition(ncid)

        end subroutine create_netcdf_parcel_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_parcel_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'pe', pe_id)

            call get_var_id(ncid, 'ke', ke_id)

            call get_var_id(ncid, 'te', te_id)

            call get_var_id(ncid, 'psi', psi_id)

            call get_var_id(ncid, 'n_parcels', npar_id)

            call get_var_id(ncid, 'n_small_parcel', nspar_id)

            call get_var_id(ncid, 'avg_lam', avg_lam_id)

            call get_var_id(ncid, 'std_lam', std_lam_id)

            call get_var_id(ncid, 'avg_vol', avg_vol_id)

            call get_var_id(ncid, 'std_vol', std_vol_id)

            call get_var_id(ncid, 'sum_vol', sum_vol_id)

            call get_var_id(ncid, 'x_rms_vorticity', rms_x_vor_id)

            call get_var_id(ncid, 'y_rms_vorticity', rms_y_vor_id)

            call get_var_id(ncid, 'z_rms_vorticity', rms_z_vor_id)

            call get_var_id(ncid, 'n_parcel_splits', n_par_split_id)

            call get_var_id(ncid, 'n_parcel_merges', n_par_merge_id)

        end subroutine read_netcdf_parcel_stats_content

        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_stats_io_timer)

            if (t <= restart_time) then
                call stop_timer(parcel_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, pe_id, pe, n_writes)
            call write_netcdf_scalar(ncid, ke_id, ke, n_writes)
            call write_netcdf_scalar(ncid, te_id, ke + pe, n_writes)
            call write_netcdf_scalar(ncid, psi_id, psi, n_writes)
            call write_netcdf_scalar(ncid, npar_id, n_parcels, n_writes)
            call write_netcdf_scalar(ncid, nspar_id, n_small, n_writes)
            call write_netcdf_scalar(ncid, avg_lam_id, avg_lam, n_writes)
            call write_netcdf_scalar(ncid, std_lam_id, std_lam, n_writes)
            call write_netcdf_scalar(ncid, avg_vol_id, avg_vol, n_writes)
            call write_netcdf_scalar(ncid, std_vol_id, std_vol, n_writes)
            call write_netcdf_scalar(ncid, sum_vol_id, sum_vol, n_writes)
            call write_netcdf_scalar(ncid, rms_x_vor_id, rms_zeta(1), n_writes)
            call write_netcdf_scalar(ncid, rms_y_vor_id, rms_zeta(2), n_writes)
            call write_netcdf_scalar(ncid, rms_z_vor_id, rms_zeta(3), n_writes)
            call write_netcdf_scalar(ncid, n_par_split_id, n_parcel_splits, n_writes)
            call write_netcdf_scalar(ncid, n_par_merge_id, n_parcel_merges, n_writes)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(parcel_stats_io_timer)

        end subroutine write_netcdf_parcel_stats
end module parcel_diagnostics_netcdf
