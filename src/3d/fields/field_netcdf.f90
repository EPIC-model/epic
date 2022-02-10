module field_netcdf
    use netcdf_utils
    use netcdf_writer
    use fields
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: netcdf_field_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: x_dim_id, y_dim_id, z_dim_id, t_dim_id

    integer            :: x_vel_id, y_vel_id, z_vel_id, &
                          x_vor_id, y_vor_id, z_vor_id, &
                          n_writes

    private :: ncid, ncfname,                           &
               x_dim_id, y_dim_id, z_dim_id, t_dim_id,  &
               x_vel_id, y_vel_id, z_vel_id,            &
               x_vor_id, y_vor_id, z_vor_id,            &
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
            call define_netcdf_dimension(ncid, "x", nx,             x_dim_id)
            call define_netcdf_dimension(ncid, "y", ny,             y_dim_id)
            call define_netcdf_dimension(ncid, "z", nz+1,           z_dim_id)
            call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, t_dim_id)

            ! define fields
            dimids = (/z_dim_id, y_dim_id, x_dim_id, t_dim_id/)

            call define_netcdf_dataset(ncid=ncid,           &
                                       name='x_velocity',   &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='m/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=x_vel_id)

            call define_netcdf_dataset(ncid=ncid,           &
                                       name='y_velocity',   &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='m/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=y_vel_id)

            call define_netcdf_dataset(ncid=ncid,           &
                                       name='z_velocity',   &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='m/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=z_vel_id)


            call define_netcdf_dataset(ncid=ncid,           &
                                       name='x_vorticity',  &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='1/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,           &
                                       name='y_vorticity',  &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='1/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,           &
                                       name='z_vorticity',  &
                                       long_name='',        &
                                       std_name='',         &
                                       unit='1/s',          &
                                       dtype=NF90_DOUBLE,   &
                                       dimids=dimids,       &
                                       varid=z_vor_id)

            call close_definition(ncid)

        end subroutine create_netcdf_field_file

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_netcdf_field_step(t, dt)
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            call start_timer(netcdf_field_timer)

! c #ifdef ENABLE_VERBOSE
!             if (verbose) then
!                 print "(a18)", "write fields to netcdf"
!             endif
! c #endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            call write_netcdf_fields

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(netcdf_field_timer)

        end subroutine write_netcdf_field_step


        ! Write field datasets (called from write_netcdf_field_step).
        ! @param[in] iter is the number of the write
        subroutine write_netcdf_fields
!             integer(hid_t)             :: group
            character(:), allocatable  :: name
            logical                    :: created
            integer                    :: cnt(4), start(4)

            ! time step to write [step(4) is the time]
            cnt   = (/ nz+1, ny, nx, 1    /)
            start = (/ 1,    1,  1,  n_writes /)

            print *, "iter", n_writes


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
!             call write_h5_dataset(ncid, name, "symmetry volume", &
!                                   sym_volg(0:nz, 0:ny-1, 0:nx-1))
!
!             call write_h5_dataset(ncid, name, "velocity gradient tensor", &
!                                   velgradg(0:nz, 0:ny-1, 0:nx-1, :))
!
!             call write_h5_dataset(ncid, name, "vorticity tendency", &
!                                   vtend(0:nz, 0:ny-1, 0:nx-1, :))
! #endif

!             call close_h5_group(group)
        end subroutine write_netcdf_fields

end module field_netcdf
