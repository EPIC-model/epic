module parcel_netcdf
    use netcdf_utils
    use netcdf_writer
    use parcel_container, only : parcels, n_parcels
    use parameters, only : nx, ny, nz, extent, lower
    implicit none

    integer :: n_writes

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id, vol_id, buo_id,  &
                          x_pos_id, y_pos_id, z_pos_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          b11_id, b12_id, b13_id,       &
                          b22_id, b23_id

#ifndef ENABLE_DRY_MODE
    integer :: hum_id
#endif

    private :: ncid, ncfname, n_writes, npar_dim_id,    &
               x_pos_id, y_pos_id, z_pos_id,            &
               x_vor_id, y_vor_id, z_vor_id,            &
               b11_id, b12_id, b13_id, b22_id, b23_id,  &
               vol_id, buo_id
#ifndef ENABLE_DRY_MODE
    private :: hum_id
#endif

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            character(:), allocatable :: name

            ncfname =  basename // '_parcels.nc'

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                ! tag the parcel file number inside
                return
            endif

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)


            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_parcels',                  &
                                         dimsize=n_parcels,                 &
                                         dimid=npar_dim_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_position',                   &
                                       long_name='x position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=(/npar_dim_id/),              &
                                       varid=x_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_position',                   &
                                       long_name='y position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=(/npar_dim_id/),              &
                                       varid=y_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_position',                   &
                                       long_name='z position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=(/npar_dim_id/),              &
                                       varid=z_pos_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B11',                              &
                                       long_name='B11 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=b11_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B12',                              &
                                       long_name='B12 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=b12_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B13',                              &
                                       long_name='B13 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=b13_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B22',                              &
                                       long_name='B22 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=b22_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B23',                              &
                                       long_name='B23 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=b23_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='volume',                           &
                                       long_name='parcel volume',               &
                                       std_name='',                             &
                                       unit='m^3',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=vol_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='x_vorticity',                      &
                                       long_name='x vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='y_vorticity',                      &
                                       long_name='y vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='z_vorticity',                      &
                                       long_name='z vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=z_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='buoyancy',                         &
                                       long_name='parcel buoyancy',             &
                                       std_name='',                             &
                                       unit='m/s^2',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=buo_id)

#ifndef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='humidity',                         &
                                       long_name='parcel humidity',             &
                                       std_name='',                             &
                                       unit='-',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=(/npar_dim_id/),                  &
                                       varid=hum_id)
#endif

            call close_definition(ncid)

        end subroutine create_netcdf_parcel_file

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcel_step(t)
            double precision, intent(in)    :: t

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

!             ! write time
!             call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            call write_netcdf_parcels

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

        end subroutine write_netcdf_parcel_step


        ! Write parcel datasets (called from write_netcdf_parcel_step).
        subroutine write_netcdf_parcels
            call write_netcdf_dataset(ncid, x_pos_id, parcels%position(1, 1:n_parcels))
            call write_netcdf_dataset(ncid, y_pos_id, parcels%position(2, 1:n_parcels))
            call write_netcdf_dataset(ncid, z_pos_id, parcels%position(3, 1:n_parcels))

            call write_netcdf_dataset(ncid, b11_id, parcels%B(1, 1:n_parcels))
            call write_netcdf_dataset(ncid, b12_id, parcels%B(2, 1:n_parcels))
            call write_netcdf_dataset(ncid, b13_id, parcels%B(3, 1:n_parcels))
            call write_netcdf_dataset(ncid, b22_id, parcels%B(4, 1:n_parcels))
            call write_netcdf_dataset(ncid, b23_id, parcels%B(5, 1:n_parcels))

            call write_netcdf_dataset(ncid, vol_id, parcels%volume(1:n_parcels))

            call write_netcdf_dataset(ncid, x_vor_id, parcels%vorticity(1, 1:n_parcels))
            call write_netcdf_dataset(ncid, y_vor_id, parcels%vorticity(2, 1:n_parcels))
            call write_netcdf_dataset(ncid, z_vor_id, parcels%vorticity(3, 1:n_parcels))

            call write_netcdf_dataset(ncid, buo_id, parcels%buoyancy(1:n_parcels))

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, hum_id, parcels%humidity(1:n_parcels))
#endif

        end subroutine write_netcdf_parcels
end module parcel_netcdf
