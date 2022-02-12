module field_netcdf
    use netcdf_utils
    use netcdf_writer
    use fields
    implicit none

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: x_dim_id, y_dim_id, z_dim_id, t_dim_id,       &
                          x_axis_id, y_axis_id, z_axis_id, t_axis_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          n_writes

    private :: ncid, ncfname,                               &
               x_dim_id, y_dim_id, z_dim_id, t_dim_id,      &
               x_axis_id, y_axis_id, z_axis_id, t_axis_id,  &
               x_vel_id, y_vel_id, z_vel_id,                &
               x_vor_id, y_vor_id, z_vor_id,                &
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
                ncerr = nf90_inquire_dimension(ncid, t_dim_id, name, n_writes)
                call check_netcdf_error("Failed to inquire the dimension.")
                call close_netcdf_file(ncid)
                return
            endif

            n_writes = 1

            call create_netcdf_file(ncfname, overwrite, ncid)

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
                dimids=(/y_dim_id/),                                       &
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
                unit='s',                                                   &
                dtype=NF90_DOUBLE,                                          &
                dimids=(/t_dim_id/),                                        &
                varid=t_axis_id)

            ncerr = nf90_put_att(ncid, t_axis_id, "axis", 't')
            call check_netcdf_error("Failed to add axis attribute.")

            ! define fields
            dimids = (/z_dim_id, y_dim_id, x_dim_id, t_dim_id/)

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

            call close_definition(ncid)

        end subroutine create_netcdf_field_file

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_step(t)
            double precision, intent(in)    :: t

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            if (n_writes == 1) then
                call write_netcdf_projected_axes
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            call write_netcdf_fields

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

        end subroutine write_netcdf_field_step


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


        ! Write field datasets (called from write_netcdf_field_step).
        subroutine write_netcdf_fields
            integer :: cnt(4), start(4)

            ! time step to write [step(4) is the time]
            cnt   = (/ nz+1, ny, nx, 1        /)
            start = (/ 1,    1,  1,  n_writes /)

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

!             call write_h5_dataset(ncid, name, "total buoyancy", &
!                                   tbuoyg(0:nz, 0:ny-1, 0:nx-1))
!
! #ifdef ENABLE_DIAGNOSE
! #ifndef ENABLE_DRY_MODE
!             call write_h5_dataset(ncid, name, "dry buoyancy", &
!                                   dbuoyg(0:nz, 0:ny-1, 0:nx-1))
! #endif
!             call write_h5_dataset(ncid, name, "volume", &
!                                   volg(0:nz, 0:ny-1, 0:nx-1))
!
!             call write_h5_dataset(ncid, name, "nparg", &
!                                   nparg(0:nz-1, :, :))
! #endif
!
! #ifndef NDEBUG
!             call write_h5_dataset(ncid, name, "velocity gradient tensor", &
!                                   velgradg(0:nz, 0:ny-1, 0:nx-1, :))
!
!             call write_h5_dataset(ncid, name, "vorticity tendency", &
!                                   vtend(0:nz, 0:ny-1, 0:nx-1, :))
! #endif

        end subroutine write_netcdf_fields

        subroutine read_netcdf_domain(fname)
            character(*), intent(in) :: fname
            integer                  :: nc_id
            double precision         :: ncells(3)

            call open_netcdf_file(ncfname, NF90_NOWRITE, nc_id)
!             call read_h5_box(nc_id, ncells, extent, lower)
            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)
            call close_netcdf_file(nc_id)

        end subroutine read_netcdf_domain

end module field_netcdf
