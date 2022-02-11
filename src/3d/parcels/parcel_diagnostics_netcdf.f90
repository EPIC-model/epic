! =============================================================================
!                      Write parcel diagnostics to NetCDF
! =============================================================================
module parcel_diagnostics_netcdf
    use netcdf_utils
    use netcdf_writer
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics
    use omp_lib
    implicit none


    private

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: t_axis_id, t_dim_id, n_writes,            &
                          pe_id, ke_id, te_id, npar_id, nspar_id,   &
                          rms_x_vor_id, rms_y_vor_id, rms_z_vor_id, &
                          avg_lam_id, std_lam_id,                   &
                          avg_vol_id, std_vol_id

    public :: create_netcdf_parcel_stats_file,  &
              write_netcdf_parcel_stats_step


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

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                ncerr = nf90_inquire_dimension(ncid, t_dim_id, name, n_writes)
                call check_netcdf_error("Failed to inquire the dimension.")
                call close_netcdf_file(ncid)
                return
            endif

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)

            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='t',                          &
                                         dimsize=NF90_UNLIMITED,            &
                                         dimid=t_dim_id)

            ncerr = nf90_put_att(ncid, t_axis_id, "axis", 't')
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='pe',                                                  &
                long_name='potential energy',                               &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=pe_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='ke',                                                  &
                long_name='kinetic energy',                                 &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=ke_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='te',                                                  &
                long_name='total energy',                                   &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=te_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='te',                                                  &
                long_name='total energy',                                   &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=te_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='n_parcels',                                           &
                long_name='number of parcels',                              &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=npar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='n_small_parcel',                                      &
                long_name='number of small parcels',                        &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=nspar_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg aspect ratio',                                    &
                long_name='average aspect ratio',                           &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_lam_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='std dev aspect ratio',                                &
                long_name='standard deviation aspect ratio',                &
                std_name='',                                                &
                unit='-',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=std_lam_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='avg volume',                                          &
                long_name='average volume',                                 &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=avg_vol_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='std dev volume',                                      &
                long_name='standard deviation volume',                      &
                std_name='',                                                &
                unit='m^3',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=std_vol_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='x rms vorticity',                                     &
                long_name='root mean square of x vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_x_vor_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='y rms vorticity',                                     &
                long_name='root mean square of y vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_y_vor_id)

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z rms vorticity',                                     &
                long_name='root mean square of z vorticity component',      &
                std_name='',                                                &
                unit='1/s',                                                 &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=rms_z_vor_id)

            call close_definition(ncid)

        end subroutine create_netcdf_parcel_stats_file

        ! Write a step in the parcel diagnostic file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_stats_step(t)
            double precision, intent(in)    :: t

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
            call write_netcdf_scalar(ncid, rms_x_vor_id, rms_zeta(1), n_writes)
            call write_netcdf_scalar(ncid, rms_y_vor_id, rms_zeta(2), n_writes)
            call write_netcdf_scalar(ncid, rms_z_vor_id, rms_zeta(3), n_writes)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

        end subroutine write_netcdf_parcel_stats_step
end module parcel_diagnostics_netcdf