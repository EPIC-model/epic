module parcel_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use parcel_container, only : parcels, n_parcels
    use parameters, only : nx, ny, nz, extent, lower, max_num_parcels
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use iomanip, only : zfill
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    integer :: parcel_io_timer

    integer :: n_writes = 1
    character(len=512) :: ncbasename

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id, vol_id, buo_id,  &
                          x_pos_id, y_pos_id, z_pos_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          b11_id, b12_id, b13_id,       &
                          b22_id, b23_id,               &
                          t_axis_id, t_dim_id
    double precision   :: restart_time

#ifndef ENABLE_DRY_MODE
    integer :: hum_id
#endif

    private :: ncid, ncfname, n_writes, npar_dim_id,    &
               x_pos_id, y_pos_id, z_pos_id,            &
               x_vor_id, y_vor_id, z_vor_id,            &
               b11_id, b12_id, b13_id, b22_id, b23_id,  &
               vol_id, buo_id, t_dim_id, t_axis_id,     &
               restart_time

#ifndef ENABLE_DRY_MODE
    private :: hum_id
#endif

    private :: ncbasename


    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: dimids(2)

            ncfname =  basename // '_' // zfill(n_writes) // '_parcels.nc'

            ncbasename = basename

            restart_time = -one

            if (l_restart) then
                ! find the last parcel file in order to set "n_writes" properly
                call exist_netcdf_file(ncfname, l_exist)
                do while (l_exist)
                    n_writes = n_writes + 1
                    ncfname =  basename // '_' // zfill(n_writes) // '_parcels.nc'
                    call exist_netcdf_file(ncfname, l_exist)
                    if (l_exist) then
                        call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                        call get_time(ncid, restart_time)
                        call close_netcdf_file(ncid)
                    endif
                enddo
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   epic_version=package_version, &
                                   file_type='parcels',          &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_parcels',                  &
                                         dimsize=n_parcels,                 &
                                         dimid=npar_dim_id)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            dimids = (/npar_dim_id, t_dim_id/)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_position',                   &
                                       long_name='x position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_position',                   &
                                       long_name='y position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_position',                   &
                                       long_name='z position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_pos_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B11',                              &
                                       long_name='B11 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b11_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B12',                              &
                                       long_name='B12 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b12_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B13',                              &
                                       long_name='B13 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b13_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B22',                              &
                                       long_name='B22 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b22_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B23',                              &
                                       long_name='B23 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b23_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='volume',                           &
                                       long_name='parcel volume',               &
                                       std_name='',                             &
                                       unit='m^3',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=vol_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='x_vorticity',                      &
                                       long_name='x vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='y_vorticity',                      &
                                       long_name='y vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='z_vorticity',                      &
                                       long_name='z vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=z_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='buoyancy',                         &
                                       long_name='parcel buoyancy',             &
                                       std_name='',                             &
                                       unit='m/s^2',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=buo_id)

#ifndef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='humidity',                         &
                                       long_name='parcel humidity',             &
                                       std_name='',                             &
                                       unit='1',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=hum_id)
#endif

            call close_definition(ncid)

        end subroutine create_netcdf_parcel_file

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcels(t)
            double precision, intent(in) :: t
            integer                      :: cnt(2), start(2)

            call start_timer(parcel_io_timer)

            if (t <= restart_time) then
                call stop_timer(parcel_io_timer)
                return
            endif

            call create_netcdf_parcel_file(trim(ncbasename), .true., .false.)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! time step to write [step(2) is the time]
            cnt   = (/ n_parcels, 1 /)
            start = (/ 1,         1 /)

            call write_netcdf_dataset(ncid, x_pos_id, parcels%position(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, y_pos_id, parcels%position(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, z_pos_id, parcels%position(3, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, b11_id, parcels%B(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b12_id, parcels%B(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b13_id, parcels%B(3, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b22_id, parcels%B(4, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b23_id, parcels%B(5, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, vol_id, parcels%volume(1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, parcels%vorticity(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, parcels%vorticity(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, parcels%vorticity(3, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, buo_id, parcels%buoyancy(1:n_parcels), start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, hum_id, parcels%humidity(1:n_parcels), start, cnt)
#endif
            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(parcel_io_timer)

        end subroutine write_netcdf_parcels

        subroutine read_netcdf_parcels(fname)
            character(*),     intent(in) :: fname
            logical                      :: l_valid = .false.
            integer                      :: cnt(2), start(2)

            call start_timer(parcel_io_timer)

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            call get_num_parcels(ncid, n_parcels)

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                stop
            endif

            ! time step to read [step(2) is the time]
            cnt   = (/ n_parcels, 1 /)
            start = (/ 1,         1 /)

            ! Be aware that the starting index of buffer_1d and buffer_2d
            ! is 0; hence, the range is 0:n_parcels-1 in contrast to the
            ! parcel container where it is 1:n_parcels.

            if (has_dataset(ncid, 'B11')) then
                call read_netcdf_dataset(ncid, 'B11', parcels%B(1, 1:n_parcels), start, cnt)
            else
                print *, "The parcel shape component B11 must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'B12')) then
                call read_netcdf_dataset(ncid, 'B12', parcels%B(2, 1:n_parcels), start, cnt)
            else
                print *, "The parcel shape component B12 must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'B13')) then
                call read_netcdf_dataset(ncid, 'B13', parcels%B(3, 1:n_parcels), start, cnt)
            else
                print *, "The parcel shape component B13 must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'B22')) then
                call read_netcdf_dataset(ncid, 'B22', parcels%B(4, 1:n_parcels), start, cnt)
            else
                print *, "The parcel shape component B22 must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'B23')) then
                call read_netcdf_dataset(ncid, 'B23', parcels%B(5, 1:n_parcels), start, cnt)
            else
                print *, "The parcel shape component B23 must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'x_position')) then
                call read_netcdf_dataset(ncid, 'x_position', &
                                         parcels%position(1, 1:n_parcels), start, cnt)
            else
                print *, "The parcel x position must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'y_position')) then
                call read_netcdf_dataset(ncid, 'y_position', &
                                         parcels%position(2, 1:n_parcels), start, cnt)
            else
                print *, "The parcel y position must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'z_position')) then
                call read_netcdf_dataset(ncid, 'z_position', &
                                         parcels%position(3, 1:n_parcels), start, cnt)
            else
                print *, "The parcel z position must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'volume')) then
                call read_netcdf_dataset(ncid, 'volume', &
                                         parcels%volume(1:n_parcels), start, cnt)
            else
                print *, "The parcel volume must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'x_vorticity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'x_vorticity', &
                                         parcels%vorticity(1, 1:n_parcels), start, cnt)
            endif

            if (has_dataset(ncid, 'y_vorticity')) then
                call read_netcdf_dataset(ncid, 'y_vorticity', &
                                         parcels%vorticity(2, 1:n_parcels), start, cnt)
            endif

            if (has_dataset(ncid, 'z_vorticity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'z_vorticity', &
                                         parcels%vorticity(3, 1:n_parcels), start, cnt)
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'buoyancy', &
                                         parcels%buoyancy(1:n_parcels), start, cnt)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'humidity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'humidity', &
                                         parcels%humidity(1:n_parcels), start, cnt)
            endif
#endif

            if (.not. l_valid) then
                print *, "Either the parcel buoyancy or vorticity must be present! Exiting."
                stop
            endif

            call close_netcdf_file(ncid)

            call stop_timer(parcel_io_timer)

        end subroutine read_netcdf_parcels

end module parcel_netcdf
