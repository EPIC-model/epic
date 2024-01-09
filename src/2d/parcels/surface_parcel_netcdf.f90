module surface_parcel_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use surface_parcel_container, only : n_top_parcels, top_parcels     &
                                       , n_bot_parcels, bot_parcels     &
                                       , surface_parcel_container_type  &
                                       , get_surface_parcel_length
    use parameters, only : nx, nz, extent, lower, max_num_parcels
    use config, only : package_version, cf_version
    use iomanip, only : zfill
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    integer :: n_writes = 1
    character(len=512) :: ncbasename

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id, len_id, buo_id,  &
                          x_pos_id, z_pos_id, vor_id,   &
                          b11_id, b12_id,               &
                          t_axis_id, t_dim_id
    double precision   :: restart_time

#ifndef ENABLE_DRY_MODE
    integer :: hum_id
#endif

    private :: ncid, ncfname, n_writes, npar_dim_id,        &
               x_pos_id, vor_id, len_id, buo_id,            &
               t_axis_id, t_dim_id,                         &
               restart_time
#ifndef ENABLE_DRY_MODE
    private :: hum_id
#endif

    private :: ncbasename, create_netcdf_parcel_file_, write_netcdf_parcels_

    contains

        ! Create the surface parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_surface_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in) :: basename
            logical,      intent(in) :: overwrite
            logical,      intent(in) :: l_restart

            call create_netcdf_parcel_file_(basename, overwrite, l_restart, n_top_parcels, 'top')
            call create_netcdf_parcel_file_(basename, overwrite, l_restart, n_bot_parcels, 'bot')

        end subroutine create_netcdf_surface_parcel_file

        ! Create the surface parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_file_(basename, overwrite, l_restart, n_par, surf)
            character(*), intent(in) :: basename
            logical,      intent(in) :: overwrite
            logical,      intent(in) :: l_restart
            integer,      intent(in) :: n_par
            character(3), intent(in) :: surf
            logical                  :: l_exist
            integer                  :: dimids(2)

            ncfname =  basename // '_' // surf // '-surface_' // zfill(n_writes) // '_parcels.nc'

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
                                   version_tag=package_version,  &
                                   file_type='parcels',          &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_parcels',                  &
                                         dimsize=n_par,                     &
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

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='length',                           &
                                       long_name='parcel length',               &
                                       std_name='',                             &
                                       unit='m',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=len_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='vorticity',                        &
                                       long_name='',                            &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=vor_id)

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

            call close_netcdf_file(ncid)

        end subroutine create_netcdf_parcel_file_

        subroutine write_netcdf_surface_parcels(t)
            double precision, intent(in) :: t

            call write_netcdf_parcels_(t, n_top_parcels, top_parcels, 'top')
            call write_netcdf_parcels_(t, n_bot_parcels, bot_parcels, 'bot')

            ! increment counter
            n_writes = n_writes + 1

        end subroutine write_netcdf_surface_parcels

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcels_(t, n_par, sp, surf)
            double precision,                    intent(in) :: t
            integer,                             intent(in) :: n_par
            type(surface_parcel_container_type), intent(in) :: sp
            character(3),                        intent(in) :: surf
            integer                                         :: cnt(2), start(2), n
            double precision                                :: length(n_par)

            if (t <= restart_time) then
                return
            endif

            call create_netcdf_parcel_file_(trim(ncbasename), .true., .false., n_par, surf)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! time step to write [step(2) is the time]
            cnt   = (/ n_par, 1 /)
            start = (/ 1,     1 /)

            do n = 1, n_par
                length(n) = get_surface_parcel_length(n, sp)
            enddo


            call write_netcdf_dataset(ncid, x_pos_id, sp%position(1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, len_id, length, start, cnt)

            call write_netcdf_dataset(ncid, vor_id, sp%vorticity(1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, buo_id, sp%buoyancy(1:n_par), start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, hum_id, sp%humidity(1:n_par), start, cnt)
#endif
            call close_netcdf_file(ncid)

        end subroutine write_netcdf_parcels_

!         subroutine read_netcdf_parcels(fname)
!             character(*),     intent(in) :: fname
!             logical                      :: l_valid = .false.
!             integer                      :: cnt(2), start(2)
!
!             call start_timer(parcel_io_timer)
!
!             call open_netcdf_file(fname, NF90_NOWRITE, ncid)
!
!             call get_num_parcels(ncid, n_parcels)
!
!             if (n_parcels > max_num_parcels) then
!                 print *, "Number of parcels exceeds limit of", &
!                           max_num_parcels, ". Exiting."
!                 stop
!             endif
!
!             ! time step to read [step(2) is the time]
!             cnt   = (/ n_parcels, 1 /)
!             start = (/ 1,         1 /)
!
!             ! Be aware that the starting index of buffer_1d and buffer_2d
!             ! is 0; hence, the range is 0:n_parcels-1 in contrast to the
!             ! parcel container where it is 1:n_parcels.
!
!             if (has_dataset(ncid, 'B11')) then
!                 call read_netcdf_dataset(ncid, 'B11', parcels%B(1, 1:n_parcels), start, cnt)
!             else
!                 print *, "The parcel shape component B11 must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'B12')) then
!                 call read_netcdf_dataset(ncid, 'B12', parcels%B(2, 1:n_parcels), start, cnt)
!             else
!                 print *, "The parcel shape component B12 must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'x_position')) then
!                 call read_netcdf_dataset(ncid, 'x_position', &
!                                          parcels%position(1, 1:n_parcels), start, cnt)
!             else
!                 print *, "The parcel x position must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'z_position')) then
!                 call read_netcdf_dataset(ncid, 'z_position', &
!                                          parcels%position(2, 1:n_parcels), start, cnt)
!             else
!                 print *, "The parcel z position must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'volume')) then
!                 call read_netcdf_dataset(ncid, 'volume', &
!                                          parcels%volume(1:n_parcels), start, cnt)
!             else
!                 print *, "The parcel volume must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'vorticity')) then
!                 l_valid = .true.
!                 call read_netcdf_dataset(ncid, 'vorticity', &
!                                          parcels%vorticity(1:n_parcels), start, cnt)
!             endif
!
!             if (has_dataset(ncid, 'buoyancy')) then
!                 l_valid = .true.
!                 call read_netcdf_dataset(ncid, 'buoyancy', &
!                                          parcels%buoyancy(1:n_parcels), start, cnt)
!             endif
!
! #ifndef ENABLE_DRY_MODE
!             if (has_dataset(ncid, 'humidity')) then
!                 l_valid = .true.
!                 call read_netcdf_dataset(ncid, 'humidity', &
!                                          parcels%humidity(1:n_parcels), start, cnt)
!             endif
! #endif
!
!             if (.not. l_valid) then
!                 print *, "Either the parcel buoyancy or vorticity must be present! Exiting."
!                 stop
!             endif
!
!             call close_netcdf_file(ncid)
!
!             call stop_timer(parcel_io_timer)
!
!         end subroutine read_netcdf_parcels

end module surface_parcel_netcdf
