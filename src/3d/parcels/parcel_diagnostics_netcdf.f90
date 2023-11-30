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
    use parcel_merging, only : n_parcel_merges
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
    integer            :: parcel_stats_io_timer

    integer, parameter :: NC_APE        = 1     &
                        , NC_KE         = 2     &
                        , NC_TE         = 3     &
                        , NC_EN         = 4     &
                        , NC_NPAR       = 5     &
                        , NC_NSPAR      = 6     &
                        , NC_RMS_X_VOR  = 7     &
                        , NC_RMS_Y_VOR  = 8     &
                        , NC_RMS_Z_VOR  = 9     &
                        , NC_AVG_LAM    = 10    &
                        , NC_STD_LAM    = 11    &
                        , NC_AVG_VOL    = 12    &
                        , NC_STD_VOL    = 13    &
                        , NC_SUM_VOL    = 14    &
                        , NC_NPAR_SPLIT = 15    &
                        , NC_NPAR_MERGE = 16    &
                        , NC_MIN_BUOY   = 17    &
                        , NC_MAX_BUOY   = 18

    type(netcdf_info) :: nc_dset(NC_MAX_BUOY)

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
            integer                   :: start(1), cnt(1), n

            if (world%rank /= world%root) then
                return
            endif

            call set_netcdf_parcel_diagnostics_info

            nc_dset(:)%l_enabled = .true.

            nc_dset(NC_APE)%l_enabled = (ape_calculation == 'ape density')

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

            ! define parcel diagnostics
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

        end subroutine create_netcdf_parcel_stats_file

        ! Pre-condition: Assumes an open file and nc_dset being initialised.
        subroutine read_netcdf_parcel_stats_content
            integer :: n

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 't', t_axis_id)

            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
                endif
            enddo

        end subroutine read_netcdf_parcel_stats_content

        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_stats(t)
            double precision, intent(in)    :: t

            call start_timer(parcel_stats_io_timer)

            if (world%rank /= world%root) then
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

            call write_diagnostic(NC_APE, parcel_stats(IDX_APE))
            call write_diagnostic(NC_KE, parcel_stats(IDX_KE))
            call write_diagnostic(NC_TE, parcel_stats(IDX_KE) + parcel_stats(IDX_APE))
            call write_diagnostic(NC_EN, parcel_stats(IDX_ENSTROPHY))
            call write_diagnostic(NC_NPAR, parcel_stats(IDX_NTOT_PAR))
            call write_diagnostic(NC_NSPAR, parcel_stats(IDX_N_SMALL))
            call write_diagnostic(NC_AVG_LAM, parcel_stats(IDX_AVG_LAM))
            call write_diagnostic(NC_STD_LAM, parcel_stats(IDX_STD_LAM))
            call write_diagnostic(NC_AVG_VOL, parcel_stats(IDX_AVG_VOL))
            call write_diagnostic(NC_STD_VOL, parcel_stats(IDX_STD_VOL))
            call write_diagnostic(NC_SUM_VOL, parcel_stats(IDX_SUM_VOL))
            call write_diagnostic(NC_RMS_X_VOR, parcel_stats(IDX_RMS_XI))
            call write_diagnostic(NC_RMS_Y_VOR, parcel_stats(IDX_RMS_ETA))
            call write_diagnostic(NC_RMS_Z_VOR, parcel_stats(IDX_RMS_ZETA))
            call write_diagnostic(NC_NPAR_SPLIT, parcel_stats(IDX_NSPLITS))
            call write_diagnostic(NC_NPAR_MERGE, parcel_stats(IDX_NMERGES))
            call write_diagnostic(NC_MIN_BUOY, parcel_stats(IDX_MIN_BUOY))
            call write_diagnostic(NC_MAX_BUOY, parcel_stats(IDX_MAX_BUOY))

            ! increment counter
            n_writes = n_writes + 1

            ! reset counters for parcel operations
            n_parcel_splits = 0
            n_parcel_merges = 0

            call close_netcdf_file(ncid, l_serial=.true.)

            call stop_timer(parcel_stats_io_timer)

        end subroutine write_netcdf_parcel_stats

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

        subroutine set_netcdf_parcel_diagnostics_info
            nc_dset(NC_APE) = netcdf_info(                                  &
                name='ape',                                                 &
                long_name='domain-averaged available potential energy',     &
                std_name='',                                                &
                unit='m^2/s^2',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_KE) = netcdf_info(                                   &
                name='ke',                                                  &
                long_name='domain-averaged kinetic energy',                 &
                std_name='',                                                &
                unit='m^2/s^2',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_TE) = netcdf_info(                                   &
                name='te',                                                  &
                long_name='domain-averaged total energy',                   &
                std_name='',                                                &
                unit='m^2/s^2',                                             &
                dtype=NF90_DOUBLE)

            nc_dset(NC_EN) = netcdf_info(                                   &
                name='en',                                                  &
                long_name='domain-averaged enstrophy',                      &
                std_name='',                                                &
                unit='1/s^2',                                               &
                dtype=NF90_DOUBLE)

            nc_dset(NC_NPAR) = netcdf_info(                                 &
                name='n_parcels',                                           &
                long_name='number of parcels',                              &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_INT)

            nc_dset(NC_NSPAR) = netcdf_info(                                &
                name='n_small_parcel',                                      &
                long_name='number of small parcels',                        &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_INT)

            nc_dset(NC_AVG_LAM) = netcdf_info(                              &
                name='avg_lam',                                             &
                long_name='average aspect ratio',                           &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE)

            nc_dset(NC_STD_LAM) = netcdf_info(                              &
                name='std_lam',                                             &
                long_name='standard deviation aspect ratio',                &
                std_name='',                                                &
                unit='1',                                                   &
                dtype=NF90_DOUBLE)

            nc_dset(NC_AVG_VOL) = netcdf_info(                              &
                name='avg_vol',                                             &
                long_name='average volume',                                 &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_STD_VOL) = netcdf_info(                              &
                name='std_vol',                                             &
                long_name='standard deviation volume',                      &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_SUM_VOL) = netcdf_info(                              &
                name='sum_vol',                                             &
                long_name='total volume',                                   &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_RMS_X_VOR) = netcdf_info(                            &
                name='x_rms_vorticity',                                     &
                long_name='root mean square of x vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_RMS_Y_VOR) = netcdf_info(                            &
                name='y_rms_vorticity',                                     &
                long_name='root mean square of y vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_RMS_Z_VOR) = netcdf_info(                            &
                name='z_rms_vorticity',                                     &
                long_name='root mean square of z vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE)

            nc_dset(NC_NPAR_SPLIT) = netcdf_info(                           &
                name='n_parcel_splits',                                     &
                 long_name='number of parcel splits since last time',       &
                 std_name='',                                               &
                 unit='1',                                                  &
                 dtype=NF90_INT)

            nc_dset(NC_NPAR_MERGE) = netcdf_info(                           &
                name='n_parcel_merges',                                     &
                 long_name='number of parcel merges since last time',       &
                 std_name='',                                               &
                 unit='1',                                                  &
                 dtype=NF90_INT)

            nc_dset(NC_MIN_BUOY) = netcdf_info(                             &
                name='min_buoyancy',                                        &
                long_name='minimum parcel buoyancy',                        &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE)

            nc_dset(NC_MAX_BUOY) = netcdf_info(                             &
                name='max_buoyancy',                                        &
                long_name='maximum parcel buoyancy',                        &
                std_name='',                                                &
                unit='m/s^2',                                               &
                dtype=NF90_DOUBLE)

        end subroutine set_netcdf_parcel_diagnostics_info

end module parcel_diagnostics_netcdf
