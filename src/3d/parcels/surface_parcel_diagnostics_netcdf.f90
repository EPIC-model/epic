! =============================================================================
!                      Write parcel diagnostics to NetCDF
! =============================================================================
module surface_parcel_diagnostics_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use surface_parcel_container, only : n_lo_surf_parcels &
                                       , n_up_surf_parcels
    use surface_parcel_diagnostics
    use parameters, only : lower, extent, nx, ny
    use config, only : package_version, cf_version
    use omp_lib
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    private

    integer :: surf_parcel_stats_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          ke_id, npar_id, nspar_id,                 &
                          rms_x_vor_id, rms_y_vor_id, rms_z_vor_id, &
                          avg_lam_id, std_lam_id,                   &
                          avg_area_id, std_area_id
    double precision   :: restart_time

    public :: create_netcdf_surface_parcel_stats_files,  &
              write_netcdf_surface_parcel_stats,         &
              surf_parcel_stats_io_timer


    contains

        subroutine create_netcdf_surface_parcel_stats_files(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

            call create_netcdf_parcel_stats_file('lo', basename, overwrite, l_restart)
            n_writes = 1
            call create_netcdf_parcel_stats_file('up', basename, overwrite, l_restart)

        end subroutine create_netcdf_surface_parcel_stats_files

        ! Create the parcel diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_stats_file(which, basename, overwrite, l_restart)
            character(2), intent(in)  :: which
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            ncfname =  basename // '_' // which // '_surf_parcel_stats.nc'

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
            call write_netcdf_info(ncid=ncid,                           &
                                   version_tag=package_version,         &
                                   file_type='surface_parcel_stats',    &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

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
                name='avg_area',                                            &
                long_name='average area',                                   &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_area_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='std_area',                                            &
                long_name='standard deviation area',                        &
                std_name='',                                                &
                unit='m^2',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=std_area_id)

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

            call close_definition(ncid)

        end subroutine create_netcdf_parcel_stats_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_parcel_stats_content

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'ke', ke_id)

            call get_var_id(ncid, 'n_parcels', npar_id)

            call get_var_id(ncid, 'n_small_parcel', nspar_id)

            call get_var_id(ncid, 'avg_lam', avg_lam_id)

            call get_var_id(ncid, 'std_lam', std_lam_id)

            call get_var_id(ncid, 'avg_area', avg_area_id)

            call get_var_id(ncid, 'std_area', std_area_id)

            call get_var_id(ncid, 'x_rms_vorticity', rms_x_vor_id)

            call get_var_id(ncid, 'y_rms_vorticity', rms_y_vor_id)

            call get_var_id(ncid, 'z_rms_vorticity', rms_z_vor_id)

        end subroutine read_netcdf_parcel_stats_content

        subroutine write_netcdf_surface_parcel_stats(t)
            double precision, intent(in) :: t

            call write_netcdf_parcel_stats(n_lo_surf_parcels, 'lo', t)
            call write_netcdf_parcel_stats(n_up_surf_parcels, 'up', t)

            ! increment counter
            n_writes = n_writes + 1

        end subroutine write_netcdf_surface_parcel_stats


        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_stats(n_par, which, t)
            integer,          intent(in) :: n_par
            character(2),     intent(in) :: which
            double precision, intent(in) :: t
            integer                      :: j

            call start_timer(surf_parcel_stats_io_timer)

            j = 1
            if (which == 'up') then
                j = 2
            endif

            if (t <= restart_time) then
                call stop_timer(surf_parcel_stats_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write diagnostics
            !
            call write_netcdf_scalar(ncid, ke_id, ke(j), n_writes)
            call write_netcdf_scalar(ncid, npar_id, n_par, n_writes)
            call write_netcdf_scalar(ncid, nspar_id, n_small(j), n_writes)
            call write_netcdf_scalar(ncid, avg_lam_id, avg_lam(j), n_writes)
            call write_netcdf_scalar(ncid, std_lam_id, std_lam(j), n_writes)
            call write_netcdf_scalar(ncid, avg_area_id, avg_area(j), n_writes)
            call write_netcdf_scalar(ncid, std_area_id, std_area(j), n_writes)
            call write_netcdf_scalar(ncid, rms_x_vor_id, rms_xi(j), n_writes)
            call write_netcdf_scalar(ncid, rms_y_vor_id, rms_eta(j), n_writes)
            call write_netcdf_scalar(ncid, rms_z_vor_id, rms_zeta(j), n_writes)


            call close_netcdf_file(ncid)

            call stop_timer(surf_parcel_stats_io_timer)

        end subroutine write_netcdf_parcel_stats
end module surface_parcel_diagnostics_netcdf
