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
    use parameters, only : lower, extent, nx, ny, nz, write_zeta_boundary_flag
    use parcel_split_mod, only : n_parcel_splits
    use parcel_merge, only : n_parcel_merges
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities, ape_calculation
    implicit none


    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          ape_id, ke_id, te_id, npar_id, nspar_id,  &
                          rms_x_vor_id, rms_y_vor_id, rms_z_vor_id, &
                          avg_lam_id, std_lam_id,                   &
                          avg_vol_id, std_vol_id, sum_vol_id,       &
                          en_id, n_par_split_id, n_par_merge_id,    &
                          min_buo_id, max_buo_id

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
            integer                   :: start(1), cnt(1)

            if (comm%rank /= comm%master) then
                return
            endif

            ncfname =  basename // '_parcel_stats.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid, l_serial=.true.)
                call get_num_steps(ncid, n_writes)
                if (n_writes > 0) then
                    call get_time(ncid, restart_time)
                    call read_netcdf_parcel_stats_content
                    start = 1
                    cnt = 1
                    call close_netcdf_file(ncid, l_serial=.true.)
                    n_writes = n_writes + 1
                    return
                else
                    call close_netcdf_file(ncid, l_serial=.true.)
                    call delete_netcdf_file(ncfname)
                endif
            endif

            call create_netcdf_file(ncfname, overwrite, ncid, l_serial=.true.)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='parcel_stats',     &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            if (ape_calculation == 'ape density') then
                call define_netcdf_dataset(                                 &
                    ncid=ncid,                                              &
                    name='ape',                                             &
                    long_name='domain-averaged available potential energy', &
                    std_name='',                                            &
                    unit='m^2/s^2',                                         &
                    dtype=NF90_DOUBLE,                                      &
                    dimids=(/t_dim_id/),                                    &
                    varid=ape_id)
            endif

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='domain-averaged kinetic energy',                 &
                std_name='',                                                &
                unit='m^2/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=ke_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='te',                                                  &
                long_name='domain-averaged total energy',                   &
                std_name='',                                                &
                unit='m^2/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=te_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='en',                                                  &
                long_name='domain-averaged enstrophy',                      &
                std_name='',                                                &
                unit='1/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=en_id)

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

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_buoyancy',                                        &
                long_name='minimum parcel buoyancy',                        &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=min_buo_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_buoyancy',                                        &
                long_name='maximum parcel buoyancy',                        &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=max_buo_id)

            call close_definition(ncid)

            call close_netcdf_file(ncid, l_serial=.true.)

        end subroutine create_netcdf_parcel_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_parcel_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            if (ape_calculation == 'ape density') then
                call get_var_id(ncid, 'ape', ape_id)
            endif

            call get_var_id(ncid, 'ke', ke_id)

            call get_var_id(ncid, 'te', te_id)

            call get_var_id(ncid, 'en', en_id)

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

            call get_var_id(ncid, 'min_buoyancy', min_buo_id)

            call get_var_id(ncid, 'max_buoyancy', max_buo_id)

        end subroutine read_netcdf_parcel_stats_content

        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_stats_io_timer)

            if (comm%rank /= comm%master) then
                return
            endif

            if (t <= restart_time) then
                call stop_timer(parcel_stats_io_timer)
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
            if (ape_calculation == 'ape density') then
                call write_netcdf_scalar(ncid, ape_id, parcel_stats(IDX_APE), n_writes, l_serial=.true.)
            endif
            call write_netcdf_scalar(ncid, ke_id, parcel_stats(IDX_KE), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, te_id, parcel_stats(IDX_KE) + parcel_stats(IDX_APE), &
                                     n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, npar_id, int(parcel_stats(IDX_NTOT_PAR)), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, nspar_id, int(parcel_stats(IDX_N_SMALL)), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, avg_lam_id, parcel_stats(IDX_AVG_LAM), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, std_lam_id, parcel_stats(IDX_STD_LAM), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, avg_vol_id, parcel_stats(IDX_AVG_VOL), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, std_vol_id, parcel_stats(IDX_STD_VOL), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, sum_vol_id, parcel_stats(IDX_SUM_VOL), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, rms_x_vor_id, parcel_stats(IDX_RMS_XI), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, rms_y_vor_id, parcel_stats(IDX_RMS_ETA), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, rms_z_vor_id, parcel_stats(IDX_RMS_ZETA), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, en_id, parcel_stats(IDX_ENSTROPHY), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, n_par_split_id, parcel_stats(IDX_NSPLITS), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, n_par_merge_id, parcel_stats(IDX_NMERGES), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, min_buo_id, parcel_stats(IDX_MIN_BUOY), n_writes, l_serial=.true.)
            call write_netcdf_scalar(ncid, max_buo_id, parcel_stats(IDX_MAX_BUOY), n_writes, l_serial=.true.)

            ! increment counter
            n_writes = n_writes + 1

            ! reset counters for parcel operations
            n_parcel_splits = 0
            n_parcel_merges = 0

            call close_netcdf_file(ncid, l_serial=.true.)

            call stop_timer(parcel_stats_io_timer)

        end subroutine write_netcdf_parcel_stats
end module parcel_diagnostics_netcdf
