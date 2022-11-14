module field_netcdf
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    implicit none

    integer :: field_io_timer

    character(len=512)  :: ncfname
    integer             :: ncid
    integer             :: dimids(3)    ! = (x, z, t)
    integer             :: coord_ids(2) ! = (x, z)
    integer             :: t_axis_id
    double precision    :: restart_time

    integer             :: x_vel_id, z_vel_id, vor_id, &
                           tbuo_id, n_writes
#ifndef ENABLE_DRY_MODE
    integer             :: dbuo_id
#endif
#ifdef ENABLE_DIAGNOSE
    integer             :: vol_id, npar_id
#endif
#ifndef NDEBUG
    integer             :: sym_vol_id
#endif

    private :: ncid, ncfname,                       &
               dimids,                              &
               coord_ids, t_axis_id,                &
               x_vel_id, z_vel_id, vor_id, tbuo_id, &
               n_writes, restart_time
#ifndef ENABLE_DRY_MODE
    private :: dbuo_id
#endif
#ifdef ENABLE_DIAGNOSE
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
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)
            call write_netcdf_box(ncid, lower, extent, (/nx, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_spatial_dimensions_2d(ncid=ncid,            &
                                                     ngps=(/nx, nz+1/),    &
                                                     dimids=dimids(1:2),   &
                                                     axids=coord_ids)

            call define_netcdf_temporal_dimension(ncid, dimids(3), t_axis_id)

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

            call get_dim_id(ncid, 'x', dimids(1))

            call get_dim_id(ncid, 'z', dimids(2))

            call get_dim_id(ncid, 't', dimids(3))

            call get_var_id(ncid, 'x', coord_ids(1))

            call get_var_id(ncid, 'z', coord_ids(2))

            call get_var_id(ncid, 't', t_axis_id)

            call get_var_id(ncid, 'x_velocity', x_vel_id)

            call get_var_id(ncid, 'z_velocity', z_vel_id)

            call get_var_id(ncid, 'vorticity', vor_id)

#ifdef ENABLE_DRY_MODE
            call get_var_id(ncid, 'buoyancy',tbuo_id)
#else
            call get_var_id(ncid, 'total_buoyancy', tbuo_id)

            call get_var_id(ncid, 'dry_buoyancy', dbuo_id)
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

            if (t <= restart_time) then
                call stop_timer(field_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_axis_2d(ncid, dimids(1:2), lower, dx, (/nx, nz/))
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


#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, dbuo_id, dbuoyg(0:nz, 0:nx-1), &
                                      start, cnt)
#endif
#ifdef ENABLE_DIAGNOSE
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

end module field_netcdf
