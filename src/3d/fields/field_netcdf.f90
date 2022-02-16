module field_netcdf
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use phys_parameters, only : glati
    implicit none

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: x_dim_id, y_dim_id, z_dim_id, t_dim_id,       &
                          x_axis_id, y_axis_id, z_axis_id, t_axis_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          tbuoy_id, dbuoy_id, lbuoy_id, &
                          n_writes

    private :: ncid, ncfname,                               &
               x_dim_id, y_dim_id, z_dim_id, t_dim_id,      &
               x_axis_id, y_axis_id, z_axis_id, t_axis_id,  &
               x_vel_id, y_vel_id, z_vel_id,                &
               x_vor_id, y_vor_id, z_vor_id,                &
               tbuoy_id, dbuoy_id, lbuoy_id,                &
               n_writes

    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: dimids(4)
            character(:), allocatable :: name

            ncfname =  basename // '_fields.nc'

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                call read_netcdf_field_content
                call close_netcdf_file(ncid)
                n_writes = n_writes + 1
                return
            endif

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_attribute(ncid=ncid, name='EPIC_version', val=package_version)
            call write_netcdf_attribute(ncid=ncid, name='file_type', val='fields')
            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))
            call write_netcdf_timestamp(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='x',                          &
                                         dimsize=nx,                        &
                                         dimid=x_dim_id)

            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='y',                          &
                                         dimsize=ny,                        &
                                         dimid=y_dim_id)

            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='z',                          &
                                         dimsize=nz+1,                      &
                                         dimid=z_dim_id)

            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='t',                          &
                                         dimsize=NF90_UNLIMITED,            &
                                         dimid=t_dim_id)


            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='x',                                                   &
                long_name='x-coordinate in projected coordinate system',    &
                std_name='projection_x_coordinate',                         &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/x_dim_id/),                                        &
                varid=x_axis_id)

            ncerr = nf90_put_att(ncid, x_axis_id, "axis", 'x')
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='y',                                                   &
                long_name='y-coordinate in projected coordinate system',    &
                std_name='projection_y_coordinate',                         &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/y_dim_id/),                                        &
                varid=y_axis_id)

            ncerr = nf90_put_att(ncid, y_axis_id, "axis", 'y')
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                         &
                ncid=ncid,                                                      &
                name='z',                                                       &
                long_name='height coordinate in projected coordinate system',   &
                std_name='height',                                              &
                unit='m',                                                       &
                dtype=NF90_DOUBLE,                                              &
                dimids=(/z_dim_id/),                                            &
                varid=z_axis_id)

            ncerr = nf90_put_att(ncid, z_axis_id, "axis", 'z')
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='t',                                                   &
                long_name='time ',                                          &
                std_name='time',                                            &
                unit='seconds since 1970-01-01',                            &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=t_axis_id)

            ncerr = nf90_put_att(ncid, t_axis_id, "axis", 't')
            call check_netcdf_error("Failed to add axis attribute.")

            ncerr = nf90_put_att(ncid, t_axis_id, "calendar", &
                                 'proleptic_gregorian')
            call check_netcdf_error("Failed to add calendear attribute.")

            ! define fields
            dimids = (/x_dim_id, y_dim_id, z_dim_id, t_dim_id/)

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
                                       unit='-',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=lbuoy_id)

            call close_definition(ncid)

        end subroutine create_netcdf_field_file

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content

            call get_dim_id(ncid, 'x', x_dim_id)

            call get_dim_id(ncid, 'y', y_dim_id)

            call get_dim_id(ncid, 'z', z_dim_id)

            call get_dim_id(ncid, 't', t_dim_id)


            call get_var_id(ncid, 'x', x_axis_id)

            call get_var_id(ncid, 'y', y_axis_id)

            call get_var_id(ncid, 'z', z_axis_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'x_velocity', x_vel_id)

            call get_var_id(ncid, 'y_velocity', y_vel_id)

            call get_var_id(ncid, 'z_velocity', z_vel_id)

            call get_var_id(ncid, 'x_vorticity', x_vor_id)

            call get_var_id(ncid, 'y_vorticity', y_vor_id)

            call get_var_id(ncid, 'z_vorticity', z_vor_id)

            call get_var_id(ncid, 'buoyancy', tbuoy_id)

            call get_var_id(ncid, 'dry_buoyancy', dbuoy_id)

            call get_var_id(ncid, 'liquid_water_content', lbuoy_id)

        end subroutine read_netcdf_field_content

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(4), start(4)

            call start_timer(field_io_timer)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_projected_axes
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

            call write_netcdf_dataset(ncid, dbuoy_id, dbuoyg(0:nz, 0:ny-1, 0:nx-1),   &
                                      start, cnt)

            call write_netcdf_dataset(ncid, lbuoy_id, glati * (tbuoyg(0:nz, 0:ny-1, 0:nx-1)     &
                                                             - dbuoyg(0:nz, 0:ny-1, 0:nx-1)),   &
                                      start, cnt)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields


        ! Write x, y, z axes (called from write_netcdf_field_step).
        subroutine write_netcdf_projected_axes
            integer          :: i
            double precision :: x_axis(0:nx-1), y_axis(0:ny-1), z_axis(0:nz)


            do i = 0, nx-1
                x_axis(i) = lower(1) + dble(i) * dx(1)
            enddo

            do i = 0, ny-1
                y_axis(i) = lower(2) + dble(i) * dx(2)
            enddo

            do i = 0, nz
                z_axis(i) = lower(3) + dble(i) * dx(3)
            enddo

            call write_netcdf_dataset(ncid, x_axis_id, x_axis)
            call write_netcdf_dataset(ncid, y_axis_id, y_axis)
            call write_netcdf_dataset(ncid, z_axis_id, z_axis)

        end subroutine write_netcdf_projected_axes

end module field_netcdf
