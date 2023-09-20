module surface_parcel_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use surface_parcel_container, only : surface_parcel_container_type      &
                                       , lo_surf_parcels, n_lo_surf_parcels &
                                       , up_surf_parcels, n_up_surf_parcels
    use parameters, only : nx, ny, nz, extent, lower, max_num_surf_parcels
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use iomanip, only : zfill
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    integer :: n_writes = 1
    character(len=512) :: ncbasename

    integer :: surf_parcel_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id, area_id, buo_id,  &
                          x_pos_id, y_pos_id,            &
                          x_vor_id, y_vor_id, z_vor_id,  &
                          b11_id, b12_id, b22_id,        &
                          t_axis_id, t_dim_id
    double precision   :: restart_time

#ifndef ENABLE_DRY_MODE
    integer :: hum_id
#endif

    private :: ncid, ncfname, n_writes, npar_dim_id,        &
               x_pos_id, y_pos_id, area_id, buo_id,         &
               x_vor_id, y_vor_id, z_vor_id,                &
               b11_id, b12_id, b22_id, t_axis_id, t_dim_id, &
               restart_time

#ifndef ENABLE_DRY_MODE
    private :: hum_id
#endif

    private :: ncbasename, create_netcdf_parcel_file, write_netcdf_parcels

    contains

        subroutine create_netcdf_surface_parcel_files(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart

            call create_netcdf_parcel_file(n_lo_surf_parcels, 'lo', &
                                           basename, overwrite, l_restart)
            n_writes = 1
            call create_netcdf_parcel_file(n_up_surf_parcels, 'up', &
                                           basename, overwrite, l_restart)

        end subroutine create_netcdf_surface_parcel_files

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_file(n_par, which, basename, overwrite, l_restart)
            integer,      intent(in)  :: n_par
            character(2), intent(in)  :: which
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: dimids(2)

            ncfname =  basename // '_' // zfill(n_writes) // '_' // which // '_surf_parcels.nc'

            ncbasename = basename

            restart_time = -one

            if (l_restart) then
                ! find the last parcel file in order to set "n_writes" properly
                call exist_netcdf_file(ncfname, l_exist)
                do while (l_exist)
                    n_writes = n_writes + 1
                    ncfname =  basename // '_' // zfill(n_writes)
                    ncfname = trim(ncfname) // '_' // which // '_surf_parcels.nc'
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
                                   file_type='surface_parcels',  &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                     &
                                         name='n_surf_parcels',         &
                                         dimsize=n_par,                 &
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
                                       name='B22',                              &
                                       long_name='B22 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b22_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='area',                             &
                                       long_name='parcel area',                 &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=area_id)

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

        subroutine write_netcdf_surface_parcels(t)
            double precision, intent(in) :: t
            call write_netcdf_parcels(lo_surf_parcels, n_lo_surf_parcels, 'lo', t)
            call write_netcdf_parcels(up_surf_parcels, n_up_surf_parcels, 'up', t)

            ! increment counter
            n_writes = n_writes + 1
        end subroutine write_netcdf_surface_parcels

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcels(s_parcels, n_par, which, t)
            type(surface_parcel_container_type), intent(in) :: s_parcels
            integer,                             intent(in) :: n_par
            character(2),                        intent(in) :: which
            double precision,                    intent(in) :: t
            integer                                         :: cnt(2), start(2)

            call start_timer(surf_parcel_io_timer)

            if (t <= restart_time) then
                call stop_timer(surf_parcel_io_timer)
                return
            endif

            call create_netcdf_parcel_file(n_par, which, trim(ncbasename), .true., .false.)

            ncfname = trim(ncbasename) // '_' // zfill(n_writes)
            ncfname = trim(ncfname) // '_' // which // '_surf_parcels.nc'

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! time step to write [step(2) is the time]
            cnt   = (/ n_par, 1 /)
            start = (/ 1,     1 /)


            call write_netcdf_dataset(ncid, x_pos_id, s_parcels%position(1, 1:n_par), start, cnt)
            call write_netcdf_dataset(ncid, y_pos_id, s_parcels%position(2, 1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, b11_id, s_parcels%B(1, 1:n_par), start, cnt)
            call write_netcdf_dataset(ncid, b12_id, s_parcels%B(2, 1:n_par), start, cnt)
            call write_netcdf_dataset(ncid, b22_id, s_parcels%B(3, 1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, area_id, s_parcels%area(1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, s_parcels%vorticity(1, 1:n_par), start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, s_parcels%vorticity(2, 1:n_par), start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, s_parcels%vorticity(3, 1:n_par), start, cnt)

            call write_netcdf_dataset(ncid, buo_id, s_parcels%buoyancy(1:n_par), start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, hum_id, s_parcels%humidity(1:n_par), start, cnt)
#endif

            call close_netcdf_file(ncid)

            call stop_timer(surf_parcel_io_timer)

        end subroutine write_netcdf_parcels

!         subroutine read_netcdf_parcels(fname)
!             character(*),     intent(in) :: fname
!             integer                      :: cnt(2), start(2)
!
!             call start_timer(surf_parcel_io_timer)
!
!             call open_netcdf_file(fname, NF90_NOWRITE, ncid)
!
!             call get_num_parcels(ncid, n_par)
!
!             if (n_par > max_num_surf_parcels) then
!                 print *, "Number of parcels exceeds limit of", &
!                           max_num_surf_parcels, ". Exiting."
!                 stop
!             endif
!
!             ! time step to read [step(2) is the time]
!             cnt   = (/ n_par, 1 /)
!             start = (/ 1,     1 /)
!
!             ! Be aware that the starting index of buffer_1d and buffer_2d
!             ! is 0; hence, the range is 0:n_par-1 in contrast to the
!             ! parcel container where it is 1:n_par.
!
!             if (has_dataset(ncid, 'B11')) then
!                 call read_netcdf_dataset(ncid, 'B11', parcels%B(1, 1:n_par), start, cnt)
!             else
!                 print *, "The parcel shape component B11 must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'B12')) then
!                 call read_netcdf_dataset(ncid, 'B12', parcels%B(2, 1:n_par), start, cnt)
!             else
!                 print *, "The parcel shape component B12 must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'B22')) then
!                 call read_netcdf_dataset(ncid, 'B22', parcels%B(3, 1:n_par), start, cnt)
!             else
!                 print *, "The parcel shape component B22 must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'x_position')) then
!                 call read_netcdf_dataset(ncid, 'x_position', &
!                                          parcels%position(1, 1:n_par), start, cnt)
!             else
!                 print *, "The parcel x position must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'y_position')) then
!                 call read_netcdf_dataset(ncid, 'y_position', &
!                                          parcels%position(2, 1:n_par), start, cnt)
!             else
!                 print *, "The parcel y position must be present! Exiting."
!                 stop
!             endif
!
!             if (has_dataset(ncid, 'area')) then
!                 call read_netcdf_dataset(ncid, 'area', &
!                                          parcels%area(1:n_par), start, cnt)
!             else
!                 print *, "The parcel area must be present! Exiting."
!                 stop
!             endif
!
!             call close_netcdf_file(ncid)
!
!             call stop_timer(surf_parcel_io_timer)
!
!         end subroutine read_netcdf_parcels

end module surface_parcel_netcdf
