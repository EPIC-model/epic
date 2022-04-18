module field_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities, glati
    implicit none

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: dimids(4)     ! = (x, y, z)
    integer            :: coord_ids(3)  ! = (x, y, z)
    integer            :: t_axis_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          tbuoy_id, n_writes
#ifndef ENABLE_DRY_MODE
    integer            :: dbuoy_id, lbuoy_id
#endif

    double precision   :: restart_time

    private :: ncid, ncfname,                   &
               dimids,                          &
               coord_ids, t_axis_id,            &
               x_vel_id, y_vel_id, z_vel_id,    &
               x_vor_id, y_vor_id, z_vor_id,    &
               tbuoy_id,                        &
               n_writes, restart_time
#ifndef ENABLE_DRY_MODE
    private :: dbuoy_id, lbuoy_id
#endif

    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            ncfname =  basename // '_fields.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call get_time(ncid, restart_time)
                call read_netcdf_field_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   epic_version=package_version, &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_spatial_dimensions_3d(ncid=ncid,                &
                                                     ncells=(/nx, ny, nz/),    &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)


            call define_netcdf_temporal_dimension(ncid, dimids(4), t_axis_id)

            ! define fields
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_velocity',                   &
                                       long_name='x velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_velocity',                   &
                                       long_name='y velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_velocity',                   &
                                       long_name='z velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vel_id)


            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_vorticity',                  &
                                       long_name='x vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_vorticity',                  &
                                       long_name='y vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_vorticity',                  &
                                       long_name='z vorticity component',   &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vor_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='total buoyancy',          &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=tbuoy_id)

#ifndef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='dry_buoyancy',                 &
                                       long_name='dry buoyancy',            &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=dbuoy_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='liquid_water_content',         &
                                       long_name='liquid-water content',    &
                                       std_name='',                         &
                                       unit='1',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=lbuoy_id)
#endif
            call close_definition(ncid)

            call close_netcdf_file(ncid)

        end subroutine create_netcdf_field_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content

            call get_dim_id(ncid, 'x', dimids(1))

            call get_dim_id(ncid, 'y', dimids(2))

            call get_dim_id(ncid, 'z', dimids(3))

            call get_dim_id(ncid, 't', dimids(4))


            call get_var_id(ncid, 'x', coord_ids(1))

            call get_var_id(ncid, 'y', coord_ids(2))

            call get_var_id(ncid, 'z', coord_ids(3))

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'x_velocity', x_vel_id)

            call get_var_id(ncid, 'y_velocity', y_vel_id)

            call get_var_id(ncid, 'z_velocity', z_vel_id)

            call get_var_id(ncid, 'x_vorticity', x_vor_id)

            call get_var_id(ncid, 'y_vorticity', y_vor_id)

            call get_var_id(ncid, 'z_vorticity', z_vor_id)

            call get_var_id(ncid, 'buoyancy', tbuoy_id)

#ifndef ENABLE_DRY_MODE
            call get_var_id(ncid, 'dry_buoyancy', dbuoy_id)

            call get_var_id(ncid, 'liquid_water_content', lbuoy_id)
#endif
        end subroutine read_netcdf_field_content

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(4), start(4)

            call start_timer(field_io_timer)

            if (t <= restart_time) then
                call stop_timer(field_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_axis_3d(ncid, dimids(1:3), lower, dx, (/nx, ny, nz/))
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            ! time step to write [step(4) is the time]
            cnt   = (/ nx, ny, nz+1, 1        /)
            start = (/ 1,  1,  1,    n_writes /)

            !
            ! write fields (do not write halo cells)
            !
            call write_netcdf_dataset(ncid, x_vel_id, velog(0:nz, 0:ny-1, 0:nx-1, 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vel_id, velog(0:nz, 0:ny-1, 0:nx-1, 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vel_id, velog(0:nz, 0:ny-1, 0:nx-1, 3), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, vortg(0:nz, 0:ny-1, 0:nx-1, 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, vortg(0:nz, 0:ny-1, 0:nx-1, 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, vortg(0:nz, 0:ny-1, 0:nx-1, 3), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, tbuoy_id, tbuoyg(0:nz, 0:ny-1, 0:nx-1),   &
                                      start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, dbuoy_id, dbuoyg(0:nz, 0:ny-1, 0:nx-1),   &
                                      start, cnt)

            call write_netcdf_dataset(ncid, lbuoy_id, glati * (tbuoyg(0:nz, 0:ny-1, 0:nx-1)     &
                                                             - dbuoyg(0:nz, 0:ny-1, 0:nx-1)),   &
                                      start, cnt)
#endif
            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields

end module field_netcdf
