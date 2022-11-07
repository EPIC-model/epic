module field_netcdf
    use constants, only : one, two
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities, glati
    use mpi_layout, only : box
    use parameters, only : write_zeta_boundary_flag
    implicit none

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: dimids(4)     ! = (x, y, z)
    integer            :: coord_ids(3)  ! = (x, y, z)
    integer            :: t_axis_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          tbuoy_id, vol_id, n_writes

#ifdef ENABLE_DIAGNOSE
    integer            :: x_vtend_id, y_vtend_id, z_vtend_id, &
                          nparg_id, nsparg_id
#endif

#ifndef ENABLE_DRY_MODE
    integer            :: dbuoy_id, lbuoy_id
#endif

    double precision   :: restart_time

    private :: ncid, ncfname,                   &
               dimids,                          &
               coord_ids, t_axis_id,            &
               x_vel_id, y_vel_id, z_vel_id,    &
               x_vor_id, y_vor_id, z_vor_id,    &
               tbuoy_id, vol_id,                &
               n_writes, restart_time

#ifdef ENABLE_DIAGNOSE
    private :: x_vtend_id, y_vtend_id, z_vtend_id, &
               nparg_id, nsparg_id
#endif

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
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_spatial_dimensions_3d(ncid=ncid,                &
                                                     ngps=(/nx, ny, nz+1/),    &
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

#ifdef ENABLE_DIAGNOSE
            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_vtend',                      &
                                       long_name='x vorticity tendency',    &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_vtend_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_vtend',                      &
                                       long_name='y vorticity tendency',    &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_vtend_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_vtend',                      &
                                       long_name='z vorticity tendency',    &
                                       std_name='',                         &
                                       unit='1/s',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_vtend_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='nparg',                            &
                                       long_name='number of parcels per cell',  &
                                       std_name='',                             &
                                       unit='1',                                &
                                       dtype=NF90_INT,                          &
                                       dimids=dimids,                           &
                                       varid=nparg_id)

            call define_netcdf_dataset(ncid=ncid,                                    &
                                       name='nsparg',                                &
                                       long_name='number of small parcels per cell', &
                                       std_name='',                                  &
                                       unit='1',                                     &
                                       dtype=NF90_INT,                               &
                                       dimids=dimids,                                &
                                       varid=nsparg_id)
#endif

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

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='volume',                       &
                                       long_name='volume',                  &
                                       std_name='',                         &
                                       unit='m^3',                          &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=vol_id)

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

#ifdef ENABLE_DIAGNOSE
            call get_var_id(ncid, 'x_vtend', x_vtend_id)

            call get_var_id(ncid, 'y_vtend', y_vtend_id)

            call get_var_id(ncid, 'z_vtend', z_vtend_id)

            call get_var_id(ncid, 'nparg', nparg_id)

            call get_var_id(ncid, 'nsparg', nsparg_id)
#endif

            call get_var_id(ncid, 'x_vorticity', x_vor_id)

            call get_var_id(ncid, 'y_vorticity', y_vor_id)

            call get_var_id(ncid, 'z_vorticity', z_vor_id)

            call get_var_id(ncid, 'buoyancy', tbuoy_id)

#ifndef ENABLE_DRY_MODE
            call get_var_id(ncid, 'dry_buoyancy', dbuoy_id)

            call get_var_id(ncid, 'liquid_water_content', lbuoy_id)
#endif
            call get_var_id(ncid, 'volume', vol_id)
        end subroutine read_netcdf_field_content

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_fields(t)
            double precision, intent(in) :: t
            integer                      :: cnt(4), start(4)
            integer                      :: lo(3), hi(3)

            call start_timer(field_io_timer)

            if (t <= restart_time) then
                call stop_timer(field_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_axis_3d(ncid, dimids(1:3), lower, dx, (/nx, ny, nz/))
                call write_zeta_boundary_flag(ncid)
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            lo = box%lo
            hi = box%hi

            ! need to add 1 since start must begin with index 1
            start(1:3) = lo + 1
            start(4) = n_writes

            cnt(1:3) = hi - lo + 1
            cnt(4)   = 1

            !
            ! write fields (do not write halo cells)
            !
            call write_netcdf_dataset(ncid, x_vel_id, velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vel_id, velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vel_id, velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3), &
                                      start, cnt)

#ifdef ENABLE_DIAGNOSE
            call write_netcdf_dataset(ncid, x_vtend_id, vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vtend_id, vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vtend_id, vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, nparg_id, nparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, nsparg_id, nparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                      start, cnt)
#endif

            call write_netcdf_dataset(ncid, x_vor_id, vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2), &
                                      start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3), &
                                      start, cnt)

            call write_netcdf_dataset(ncid, tbuoy_id, tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),   &
                                      start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, dbuoy_id, dbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),   &
                                      start, cnt)

            call write_netcdf_dataset(ncid, lbuoy_id, glati * (tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))    &
                                                             - dbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))),  &
                                      start, cnt)
#endif

            call write_netcdf_dataset(ncid, vol_id, volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                      start, cnt)

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields

end module field_netcdf
