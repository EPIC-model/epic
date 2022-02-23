! =============================================================================
!                      Write parcel diagnostics to NetCDF
! =============================================================================
module parcel_diagnostics_netcdf
    use constants, only : one
    use parcel_diagnostics
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics
    use parameters, only : lower, extent, nx, nz
    use config, only : package_version
    use omp_lib
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    implicit none

    private

    integer :: parcel_stats_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          pe_id, ke_id, te_id, npar_id, nspar_id,   &
                          rms_vor_id,                               &
                          avg_lam_id, std_lam_id,                   &
                          avg_vol_id, std_vol_id,                   &
                          min_buo_id, max_buo_id,                   &
                          min_vor_id, max_vor_id
    double precision   :: restart_time

#ifdef ENABLE_DIAGNOSE
    integer             :: xb_bar_id, x2b_bar_id, xzb_bar_id,       &
                           zb_bar_id, z2b_bar_id,                   &
                           xv_bar_id, x2v_bar_id, xzv_bar_id,       &
                           zv_bar_id, z2v_bar_id
#endif

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
            character(:), allocatable :: name

            ncfname =  basename // '_parcel_stats.nc'

            call exist_netcdf_file(ncfname, l_exist)

            restart_time = -one
            n_writes = 1

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
            call write_netcdf_attribute(ncid=ncid, name='EPIC_version', val=package_version)
            call write_netcdf_attribute(ncid=ncid, name='file_type', val='parcel_stats')
            call write_netcdf_attribute(ncid=ncid, name='Conventions', val='CF-1.8')
            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))
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
                name='pe',                                                  &
                long_name='potential energy',                               &
                std_name='',                                                &
                unit='m^4/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=pe_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='kinetic energy',                                 &
                std_name='',                                                &
                unit='m^4/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=ke_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='te',                                                  &
                long_name='total energy',                                   &
                std_name='',                                                &
                unit='m^4/s^2',                                             &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=te_id)

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
                name='rms_vorticity',                                       &
                long_name='root mean square of vorticity',                  &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_vor_id)

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

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='min_vorticity',                                       &
                long_name='minimum parcel vorticity',                       &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=min_vor_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='max_vorticity',                                       &
                long_name='maximum parcel vorticity',                       &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=max_vor_id)

#ifdef ENABLE_DIAGNOSE
            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='xb_bar',                                              &
                long_name='buoyancy-weighted centre of mass in x',          &
                std_name='',                                                &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=xb_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='x2b_bar',                                             &
                long_name='buoyancy-weighted variance in x',                &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=x2b_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='zb_bar',                                              &
                long_name='buoyancy-weighted centre of mass in z',          &
                std_name='',                                                &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=zb_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z2b_bar',                                             &
                long_name='buoyancy-weighted variance in z',                &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=z2b_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='xzb_bar',                                             &
                long_name='buoyancy-weighted correlation in xz,             &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=xzb_bar_id)


            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='xv_bar',                                              &
                long_name='vorticity-weighted centre of mass in x',         &
                std_name='',                                                &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=xv_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='x2v_bar',                                             &
                long_name='vorticity-weighted variance in x',               &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=x2v_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='zv_bar',                                              &
                long_name='vorticity-weighted centre of mass in z',         &
                std_name='',                                                &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=zv_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z2v_bar',                                             &
                long_name='vorticity-weighted variance in z',               &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=z2v_bar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='xzv_bar',                                             &
                long_name='vorticity-weighted correlation in xz,            &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=xzv_bar_id)
#endif

            call close_definition(ncid)

        end subroutine create_netcdf_parcel_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_parcel_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'pe', pe_id)

            call get_var_id(ncid, 'ke', ke_id)

            call get_var_id(ncid, 'te', te_id)

            call get_var_id(ncid, 'n_parcels', npar_id)

            call get_var_id(ncid, 'n_small_parcel', nspar_id)

            call get_var_id(ncid, 'avg_lam', avg_lam_id)

            call get_var_id(ncid, 'std_lam', std_lam_id)

            call get_var_id(ncid, 'avg_vol', avg_vol_id)

            call get_var_id(ncid, 'std_vol', std_vol_id)

            call get_var_id(ncid, 'rms_vorticity', rms_vor_id)

            call get_var_id(ncid, 'min_buoyancy', min_buo_id)

            call get_var_id(ncid, 'max_buoyancy', max_buo_id)

            call get_var_id(ncid, 'min_vorticity', min_vor_id)

            call get_var_id(ncid, 'max_vorticity', max_vor_id)

#ifdef ENABLE_DIAGNOSE
            call get_var_id(ncid, 'xb_bar', xb_bar_id)

            call get_var_id(ncid, 'x2b_bar', x2b_bar_id)

            call get_var_id(ncid, 'zb_bar', zb_bar_id)

            call get_var_id(ncid, 'z2b_bar', z2b_bar_id)

            call get_var_id(ncid, 'xzb_bar', xzb_bar_id)

            call get_var_id(ncid, 'xv_bar', xv_bar_id)

            call get_var_id(ncid, 'x2v_bar', x2v_bar_id)

            call get_var_id(ncid, 'zv_bar', zv_bar_id)

            call get_var_id(ncid, 'z2v_bar', z2v_bar_id)

            call get_var_id(ncid, 'xzv_bar', xzv_bar_id)
#endif

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
            call write_netcdf_scalar(ncid, npar_id, n_parcels, n_writes)
            call write_netcdf_scalar(ncid, nspar_id, n_small, n_writes)
            call write_netcdf_scalar(ncid, avg_lam_id, avg_lam, n_writes)
            call write_netcdf_scalar(ncid, std_lam_id, std_lam, n_writes)
            call write_netcdf_scalar(ncid, avg_vol_id, avg_vol, n_writes)
            call write_netcdf_scalar(ncid, std_vol_id, std_vol, n_writes)
            call write_netcdf_scalar(ncid, rms_vor_id, rms_zeta, n_writes)
            call write_netcdf_scalar(ncid, min_buo_id, bmin, n_writes)
            call write_netcdf_scalar(ncid, max_buo_id, bmax, n_writes)
            call write_netcdf_scalar(ncid, min_vor_id, vormin, n_writes)
            call write_netcdf_scalar(ncid, max_vor_id, vormax, n_writes)

#ifdef ENABLE_DIAGNOSE
            call write_netcdf_scalar(ncid, xb_bar_id, xb_bar, n_writes)
            call write_netcdf_scalar(ncid, x2b_bar_id, x2b_bar, n_writes)
            call write_netcdf_scalar(ncid, zb_bar_id, zb_bar, n_writes)
            call write_netcdf_scalar(ncid, z2b_bar_id, z2b_bar, n_writes)
            call write_netcdf_scalar(ncid, xzb_bar_id, xzb_bar, n_writes)

            call write_netcdf_scalar(ncid, xv_bar_id, xv_bar, n_writes)
            call write_netcdf_scalar(ncid, x2v_bar_id, x2v_bar, n_writes)
            call write_netcdf_scalar(ncid, zv_bar_id, zv_bar, n_writes)
            call write_netcdf_scalar(ncid, z2v_bar_id, z2v_bar, n_writes)
            call write_netcdf_scalar(ncid, xzv_bar_id, xzv_bar, n_writes)
#endif

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(parcel_stats_io_timer)

        end subroutine write_netcdf_parcel_stats
end module parcel_diagnostics_netcdf
