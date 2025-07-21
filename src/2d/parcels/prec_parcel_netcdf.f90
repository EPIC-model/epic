module prec_parcel_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use precipitation_parcels, only : prec_parcels, n_prec_parcels
    use parameters, only : nx, nz, extent, lower, max_num_prec_parcels
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use iomanip, only : zfill
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    integer :: n_writes = 1
    character(len=512) :: ncbasename

    integer :: prec_parcel_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id, vol_id, qr_id,  &
                          x_pos_id, z_pos_id, Nr_id,   &
                          t_axis_id, t_dim_id
    double precision   :: restart_time

    private :: ncid, ncfname, n_writes, npar_dim_id,        &
               x_pos_id, z_pos_id, vol_id, qr_id, Nr_id,  &
               t_axis_id, t_dim_id,         &
               restart_time

    private :: ncbasename

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_prec_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: dimids(2)

            ncfname =  basename // '_' // zfill(n_writes) // '_prec_parcels.nc'

            ncbasename = basename

            restart_time = -one

            if (l_restart) then
                ! find the last parcel file in order to set "n_writes" properly
                call exist_netcdf_file(ncfname, l_exist)
                do while (l_exist)
                    n_writes = n_writes + 1
                    ncfname =  basename // '_' // zfill(n_writes) // '_prec_parcels.nc'
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
                                   version_tag=package_version,  &
                                   file_type='prec_parcels',          &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_prec_parcels',                  &
                                         dimsize=n_prec_parcels,                 &
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
                                       name='z_position',                   &
                                       long_name='z position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_pos_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='volume',                           &
                                       long_name='parcel volume',               &
                                       std_name='',                             &
                                       unit='m^3',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=vol_id)


            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='qr',                         &
                                       long_name='parcel rain mixing ratio',             &
                                       std_name='',                             &
                                       unit='kg/kg',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=qr_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='Nr',                               &
                                       long_name='parcel rain number',          &
                                       std_name='',                             &
                                       unit='/kg',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=Nr_id)

            call close_definition(ncid)

        end subroutine create_netcdf_prec_parcel_file

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_prec_parcels(t)
            double precision, intent(in) :: t
            integer                      :: cnt(2), start(2)

            call start_timer(prec_parcel_io_timer)

            if (t <= restart_time) then
                call stop_timer(prec_parcel_io_timer)
                return
            endif

            call create_netcdf_prec_parcel_file(trim(ncbasename), .true., .false.)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! time step to write [step(2) is the time]
            cnt   = (/ n_prec_parcels, 1 /)
            start = (/ 1,         1 /)


            call write_netcdf_dataset(ncid, x_pos_id, prec_parcels%position(1, 1:n_prec_parcels), start, cnt)
            call write_netcdf_dataset(ncid, z_pos_id, prec_parcels%position(2, 1:n_prec_parcels), start, cnt)

            call write_netcdf_dataset(ncid, vol_id, prec_parcels%volume(1:n_prec_parcels), start, cnt)

            call write_netcdf_dataset(ncid, qr_id, prec_parcels%qr(1:n_prec_parcels), start, cnt)

            call write_netcdf_dataset(ncid, Nr_id, prec_parcels%Nr(1:n_prec_parcels), start, cnt)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(prec_parcel_io_timer)

        end subroutine write_netcdf_prec_parcels

        subroutine read_netcdf_prec_parcels(fname)
            character(*),     intent(in) :: fname
            integer                      :: cnt(2), start(2)

            call start_timer(prec_parcel_io_timer)

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            call get_num_parcels(ncid, n_prec_parcels)

            if (n_prec_parcels > max_num_prec_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_prec_parcels, ". Exiting."
                stop
            endif

            ! time step to read [step(2) is the time]
            cnt   = (/ n_prec_parcels, 1 /)
            start = (/ 1,         1 /)

            ! Be aware that the starting index of buffer_1d and buffer_2d
            ! is 0; hence, the range is 0:n_prec_parcels-1 in contrast to the
            ! parcel container where it is 1:n_prec_parcels.

            if (has_dataset(ncid, 'x_position')) then
                call read_netcdf_dataset(ncid, 'x_position', &
                                         prec_parcels%position(1, 1:n_prec_parcels), start, cnt)
            else
                print *, "The parcel x position must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'z_position')) then
                call read_netcdf_dataset(ncid, 'z_position', &
                                         prec_parcels%position(2, 1:n_prec_parcels), start, cnt)
            else
                print *, "The parcel z position must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'volume')) then
                call read_netcdf_dataset(ncid, 'volume', &
                                         prec_parcels%volume(1:n_prec_parcels), start, cnt)
            else
                print *, "The parcel volume must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'qr')) then
                call read_netcdf_dataset(ncid, 'qr', &
                                         prec_parcels%qr(1:n_prec_parcels), start, cnt)
            else
                print *, "The parcel qr must be present! Exiting."
                stop
            endif

            if (has_dataset(ncid, 'Nr')) then
                call read_netcdf_dataset(ncid, 'Nr', &
                                         prec_parcels%Nr(1:n_prec_parcels), start, cnt)
            else
                print *, "The parcel Nr must be present! Exiting."
                stop
            endif

            call close_netcdf_file(ncid)

            call stop_timer(prec_parcel_io_timer)

        end subroutine read_netcdf_prec_parcels

end module prec_parcel_netcdf
