module field_netcdf
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    implicit none

    integer :: field_io_timer

    character(len=512)  :: ncfname
    integer             :: ncid
    integer             :: x_dim_id, z_dim_id, t_dim_id,       &
                           x_axis_id, z_axis_id, t_axis_id

    integer             :: x_vel_id, z_vel_id, vor_id, &
                           tbuo_id, n_writes
#ifdef ENABLE_DIAGNOSE
#ifndef ENABLE_DRY_MODE
    integer             :: dbuo_id
#endif
    integer             :: vol_id, npar_id
#endif
#ifndef NDEBUG
    integer             :: sym_vol_id
#endif

    private :: ncid, ncfname,                       &
               x_dim_id, z_dim_id, t_dim_id,        &
               x_axis_id, z_axis_id, t_axis_id,     &
               x_vel_id, z_vel_id, vor_id, tbuo_id, &
               n_writes
#ifdef ENABLE_DIAGNOSE
#ifndef ENABLE_DRY_MODE
    private :: dbuo_id
#endif
    private :: vol_id, npar_id
#endif
#ifndef NDEBUG
    private :: sym_vol_id
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
            integer                   :: dimids(3)
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
            call write_netcdf_attribute(ncid=ncid, name='Conventions', val='CF-1.9')
            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))
            call write_netcdf_timestamp(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='x',                          &
                                         dimsize=nx,                        &
                                         dimid=x_dim_id)

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
                long_name='width coordinate',                               &
                std_name='width',                                           &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/x_dim_id/),                                        &
                varid=x_axis_id)

            ncerr = nf90_put_att(ncid, x_axis_id, "axis", 'X')
            call check_netcdf_error("Failed to add axis attribute.")

            call define_netcdf_dataset(                                     &
                ncid=ncid,                                                  &
                name='z',                                                   &
                long_name='height coordinate',                              &
                std_name='height',                                          &
                unit='m',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/z_dim_id/),                                        &
                varid=z_axis_id)

            ncerr = nf90_put_att(ncid, z_axis_id, "axis", 'Z')
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

            ncerr = nf90_put_att(ncid, t_axis_id, "axis", 'T')
            call check_netcdf_error("Failed to add axis attribute.")

            ncerr = nf90_put_att(ncid, t_axis_id, "calendar", &
                                 'proleptic_gregorian')
            call check_netcdf_error("Failed to add calendear attribute.")

            ! define fields
            dimids = (/x_dim_id, z_dim_id, t_dim_id/)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_velocity',                   &
                                       long_name='x velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vel_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_velocity',                   &
                                       long_name='z velocity component',    &
                                       std_name='',                         &
                                       unit='m/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vel_id)


            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='vorticity',                    &
                                       long_name='vorticity',               &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=vor_id)

#ifdef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='buoyancy',                &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=tbuo_id)
#else
#ifdef ENABLE_DIAGNOSE
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='total_buoyancy',               &
                                       long_name='total buoyancy',          &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=tbuo_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='dry_buoyancy',                 &
                                       long_name='dry buoyancy',            &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=dbuo_id)
#endif
#endif

#ifdef ENABLE_DIAGNOSE
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='volume',                       &
                                       long_name='volume',                  &
                                       std_name='',                         &
                                       unit='m^2',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=vol_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='npar',                             &
                                       long_name='number of parcels per cell',  &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=npar_id)
#endif

#ifndef NDEBUG
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='sym_vol',                      &
                                       long_name='symmetry volume',         &
                                       std_name='',                         &
                                       unit='m^2',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=sym_vol_id)
#endif

            call close_definition(ncid)

        end subroutine create_netcdf_field_file


        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content

            call get_dim_id(ncid, 'x', x_dim_id)

            call get_dim_id(ncid, 'z', z_dim_id)

            call get_dim_id(ncid, 't', t_dim_id)

            call get_var_id(ncid, 'x', x_axis_id)

            call get_var_id(ncid, 'z', z_axis_id)

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'x_velocity', x_vel_id)

            call get_var_id(ncid, 'z_velocity', z_vel_id)

            call get_var_id(ncid, 'vorticity', vor_id)

#ifdef ENABLE_DRY_MODE
            call get_var_id(ncid, 'buoyancy',tbuo_id)
#else
#ifdef ENABLE_DIAGNOSE
            call get_var_id(ncid, 'total_buoyancy', tbuo_id)

            call get_var_id(ncid, 'dry_buoyancy', dbuo_id)
#endif
#endif

#ifdef ENABLE_DIAGNOSE
            call get_var_id(ncid, 'volume', vol_id)

            call get_var_id(ncid, 'npar', npar_id)
#endif

#ifndef NDEBUG
            call get_var_id(ncid, 'sym_vol', sym_vol_id)
#endif

        end subroutine read_netcdf_field_content

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(3), start(3)

            call start_timer(field_io_timer)

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_projected_axes
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            ! time step to write [step(4) is the time]
            cnt   = (/ nx, nz+1, 1        /)
            start = (/ 1,  1,    n_writes /)

            !
            ! write fields (do not write halo cells)
            !
            call write_netcdf_dataset(ncid, x_vel_id, velog(0:nz, 0:nx-1, 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vel_id, velog(0:nz, 0:nx-1, 2), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, vor_id, vortg(0:nz, 0:nx-1), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, tbuo_id, tbuoyg(0:nz, 0:nx-1), &
                                      start, cnt)


#ifdef ENABLE_DIAGNOSE
#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, dbuo_id, dbuoyg(0:nz, 0:nx-1), &
                                      start, cnt)
#endif
            call write_netcdf_dataset(ncid, vol_id, volg(0:nz, 0:nx-1))

            call write_netcdf_dataset(ncid, npar_id, nparg(0:nz-1, :))
#endif

#ifndef NDEBUG
            call write_netcdf_dataset(ncid, sym_vol_id, sym_volg(0:nz, 0:nx-1))
#endif


            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields


        ! Write x, y, z axes (called from write_netcdf_field_step).
        subroutine write_netcdf_projected_axes
            integer          :: i
            double precision :: x_axis(0:nx-1), z_axis(0:nz)


            do i = 0, nx-1
                x_axis(i) = lower(1) + dble(i) * dx(1)
            enddo

            do i = 0, nz
                z_axis(i) = lower(2) + dble(i) * dx(2)
            enddo

            call write_netcdf_dataset(ncid, x_axis_id, x_axis)
            call write_netcdf_dataset(ncid, z_axis_id, z_axis)

        end subroutine write_netcdf_projected_axes

end module field_netcdf
