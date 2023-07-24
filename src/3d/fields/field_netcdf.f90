module field_netcdf
    use constants, only : one, two
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use fields
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use options, only : write_netcdf_options, output
    use physics, only : write_physical_quantities, glati
    use mpi_layout, only : box
    use parameters, only : write_zeta_boundary_flag
    use mpi_utils, only : mpi_stop
    implicit none

    private

    integer :: field_io_timer

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: dimids(4)     ! = (x, y, z)
    integer            :: coord_ids(3)  ! = (x, y, z)
    integer            :: t_axis_id


    type netcdf_field_info
        character(64)  :: name      = ''
        character(128) :: long_name = ''
        character(128) :: std_name  = ''
        character(16)  :: unit      = ''
        integer        :: dtype     = -1
        integer        :: varid     = -1
        logical        :: l_enabled = .false.
    end type netcdf_field_info

    integer, parameter :: NC_X_VEL   = 1        &
                        , NC_Y_VEL   = 2        &
                        , NC_Z_VEL   = 3        &
                        , NC_X_VOR   = 4        &
                        , NC_Y_VOR   = 5        &
                        , NC_Z_VOR   = 6        &
                        , NC_X_VTEND = 7        &
                        , NC_Y_VTEND = 8        &
                        , NC_Z_VTEND = 9        &
                        , NC_NPARG   = 10       &
                        , NC_NSPARG  = 11       &
                        , NC_TBUOY   = 12       &
                        , NC_VOL     = 13

#ifndef ENABLE_DRY_MODE
    integer, parameter :: NC_DBUOY   = 14       &
                        , NC_HUM     = 15       &
                        , NC_LBUOY   = 16
#endif

    integer :: n_writes

    type(netcdf_field_info) :: nc_dset(16)

    double precision   :: restart_time

    public :: create_netcdf_field_file  &
            , write_netcdf_fields       &
            , read_netcdf_fields        &
            , field_io_timer

    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: n

            call set_netcdf_field_output

            ncfname =  basename // '_fields.nc'

            restart_time = -one
            n_writes = 1

            call exist_netcdf_file(ncfname, l_exist)

            if (l_restart .and. l_exist) then
                call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                call get_num_steps(ncid, n_writes)
                if (n_writes > 0) then
                    call get_time(ncid, restart_time)
                    call read_netcdf_field_content
                    call close_netcdf_file(ncid)
                    n_writes = n_writes + 1
                    return
                else
                    call close_netcdf_file(ncid)
                    if (world%rank == world%root) then
                        call delete_netcdf_file(ncfname)
                    endif
                endif
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
            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call define_netcdf_dataset(ncid=ncid,                       &
                                               name=nc_dset(n)%name,            &
                                               long_name=nc_dset(n)%long_name,  &
                                               std_name=nc_dset(n)%std_name,    &
                                               unit=nc_dset(n)%unit,            &
                                               dtype=nc_dset(n)%dtype,          &
                                               dimids=dimids,                   &
                                               varid=nc_dset(n)%varid)

                endif
            enddo

            call close_definition(ncid)

            call close_netcdf_file(ncid)

        end subroutine create_netcdf_field_file

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Pre-condition: Assumes an open file
        subroutine read_netcdf_field_content
            integer :: n

            call get_dim_id(ncid, 'x', dimids(1))

            call get_dim_id(ncid, 'y', dimids(2))

            call get_dim_id(ncid, 'z', dimids(3))

            call get_dim_id(ncid, 't', dimids(4))

            call get_var_id(ncid, 'x', coord_ids(1))

            call get_var_id(ncid, 'y', coord_ids(2))

            call get_var_id(ncid, 'z', coord_ids(3))

            call get_var_id(ncid, 't', t_axis_id)

            do n = 1, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
                endif
            enddo

        end subroutine read_netcdf_field_content

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

            lo = box%lo
            hi = box%hi

            ! need to add 1 since start must begin with index 1
            start(1:3) = lo + 1
            start(4) = n_writes

            cnt(1:3) = hi - lo + 1
            cnt(4)   = 1

            if (n_writes == 1) then
                call write_netcdf_axis_3d(ncid, dimids(1:3), box%lower, dx, &
                                          box%size, start(1:3), cnt(1:3))
                call write_zeta_boundary_flag(ncid)
            endif

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, n_writes)

            !
            ! write fields (do not write halo cells)
            !
            if (nc_dset(NC_X_VEL)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_X_VEL)%varid,                    &
                                          velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Y_VEL)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Y_VEL)%varid,                    &
                                          velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Z_VEL)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Z_VEL)%varid,                    &
                                          velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_X_VTEND)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_X_VTEND)%varid,                  &
                                          vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Y_VTEND)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Y_VTEND)%varid,                  &
                                          vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Z_VTEND)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Z_VTEND)%varid,                  &
                                          vtend(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_NPARG)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_NPARG)%varid,                    &
                                          nparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),     &
                                          start, cnt)
            endif

            if (nc_dset(NC_NSPARG)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_NSPARG)%varid,                   &
                                          nsparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),    &
                                          start, cnt)
            endif

            if (nc_dset(NC_X_VOR)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_X_VOR)%varid,                    &
                                          vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Y_VOR)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Y_VOR)%varid,                    &
                                          vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_Z_VOR)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_Z_VOR)%varid,                    &
                                          vortg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3),  &
                                          start, cnt)
            endif


            if (nc_dset(NC_TBUOY)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_TBUOY)%varid,                    &
                                          tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),    &
                                          start, cnt)
            endif

#ifndef ENABLE_DRY_MODE
            if (nc_dset(NC_DBUOY)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_DBUOY)%varid,                    &
                                          dbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),    &
                                         start, cnt)
            endif

            if (nc_dset(NC_LBUOY)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_LBUOY)%varid,                            &
                                          glati * (tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))    &
                                                 - dbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))),  &
                                          start, cnt)
            endif

            if (nc_dset(NC_HUM)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_HUM)%varid,                   &
                                          humg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),   &
                                          start, cnt)
            endif
#endif
            if (nc_dset(NC_VOL)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_VOL)%varid,                   &
                                          volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)),   &
                                          start, cnt)
            endif

            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(field_io_timer)

        end subroutine write_netcdf_fields

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine read_netcdf_fields(fname, step)
            character(*), intent(in) :: fname
            integer,      intent(in) :: step
            integer                  :: n_steps, lid, start(4), cnt(4), st

            call open_netcdf_file(fname, NF90_NOWRITE, lid)

            call get_num_steps(lid, n_steps)

            st = step

            if (st == -1) then
                st = n_steps
            else if ((st == 0) .or. (st > n_steps)) then
                call mpi_stop("Step number in NetCDF field file is out of bounds.")
            endif

            start(1:3) = box%lo + 1      ! need to add 1 since start must begin with index 1
            cnt(1:3) = box%size
            cnt(4) = 1
            start(4) = st

            if (has_dataset(lid, nc_dset(NC_X_VOR)%name)) then
                call read_netcdf_dataset(lid,                       &
                                         nc_dset(NC_X_VOR)%name,    &
                                         vortg(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1), &
                                               1),                  &
                                         start,                     &
                                         cnt)
            endif

            if (has_dataset(lid, nc_dset(NC_Y_VOR)%name)) then
                call read_netcdf_dataset(lid,                       &
                                         nc_dset(NC_Y_VOR)%name,    &
                                         vortg(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1), &
                                               2),                  &
                                         start,                     &
                                         cnt)
            endif

            if (has_dataset(lid, nc_dset(NC_Z_VOR)%name)) then
                call read_netcdf_dataset(lid,                       &
                                         nc_dset(NC_Z_VOR)%name,    &
                                         vortg(box%lo(3):box%hi(3), &
                                               box%lo(2):box%hi(2), &
                                               box%lo(1):box%hi(1), &
                                               3),                  &
                                         start,                     &
                                         cnt)
            endif

            if (has_dataset(lid, nc_dset(NC_TBUOY)%name)) then
                call read_netcdf_dataset(lid,                         &
                                         nc_dset(NC_TBUOY)%name,      &
                                         tbuoyg(box%lo(3):box%hi(3),  &
                                                box%lo(2):box%hi(2),  &
                                                box%lo(1):box%hi(1)), &
                                         start,                       &
                                         cnt)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(lid, nc_dset(NC_HUM)%name)) then
                call read_netcdf_dataset(lid,                       &
                                         nc_dset(NC_HUM)%name,      &
                                         humg(box%lo(3):box%hi(3),  &
                                              box%lo(2):box%hi(2),  &
                                              box%lo(1):box%hi(1)), &
                                         start,                     &
                                         cnt)
            endif
#endif
            call close_netcdf_file(lid)

        end subroutine read_netcdf_fields

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_output
            integer :: n

            call set_netcdf_field_info

            nc_dset(:)%l_enabled = .false.

            do n = 1, size(output%field_list)
                select case (output%field_list(n))
                    case ('x-vorticity')
                        nc_dset(NC_X_VOR)%l_enabled = .true.
                    case ('y-vorticity')
                        nc_dset(NC_Y_VOR)%l_enabled = .true.
                    case ('z-vorticity')
                        nc_dset(NC_Z_VOR)%l_enabled = .true.
                    case ('x-velocity')
                        nc_dset(NC_X_VEL)%l_enabled = .true.
                    case ('y-velocity')
                        nc_dset(NC_Y_VEL)%l_enabled = .true.
                    case ('z-velocity')
                        nc_dset(NC_Z_VEL)%l_enabled = .true.
                    case ('x-vorticity-tendency')
                        nc_dset(NC_X_VTEND)%l_enabled = .true.
                    case ('y-vorticity-tendency')
                        nc_dset(NC_Y_VTEND)%l_enabled = .true.
                    case ('z-vorticity-tendency')
                        nc_dset(NC_Z_VTEND)%l_enabled = .true.
                    case ('nparg')
                        nc_dset(NC_NPARG)%l_enabled = .true.
                    case ('nsparg')
                        nc_dset(NC_NSPARG)%l_enabled = .true.
                    case ('buoyancy')
                        nc_dset(NC_TBUOY)%l_enabled = .true.
#ifndef ENABLE_DRY_MODE
                    case ('dry-buoyancy')
                        nc_dset(NC_DBUOY)%l_enabled = .true.
                    case ('humidity')
                        nc_dset(NC_HUM)%l_enabled = .true.
                    case ('liquid-water-content')
                        nc_dset(NC_LBUOY)%l_enabled = .true.
#endif
                    case ('volume')
                        nc_dset(NC_VOL)%l_enabled = .true.
                    case ('all')
                        nc_dset(:)%l_enabled = .true.
                    case ('default')
                        nc_dset(NC_X_VOR)%l_enabled = .true.
                        nc_dset(NC_Y_VOR)%l_enabled = .true.
                        nc_dset(NC_Z_VOR)%l_enabled = .true.
                        nc_dset(NC_X_VEL)%l_enabled = .true.
                        nc_dset(NC_Y_VEL)%l_enabled = .true.
                        nc_dset(NC_Z_VEL)%l_enabled = .true.
                        nc_dset(NC_TBUOY)%l_enabled = .true.
#ifndef ENABLE_DRY_MODE
                        nc_dset(NC_DBUOY)%l_enabled = .true.
                        nc_dset(NC_HUM)%l_enabled   = .true.
                        nc_dset(NC_LBUOY)%l_enabled = .true.
#endif
                        nc_dset(NC_VOL)%l_enabled   = .true.
                    case ('')
                        ! do nothing
                    case default
                        call mpi_stop("No field '" // output%field_list(n) // "' available.")
                end select
            enddo

        end subroutine set_netcdf_field_output

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_field_info
            nc_dset(NC_X_VEL) = netcdf_field_info(name='x_velocity',                    &
                                                  long_name='x velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VEL) = netcdf_field_info(name='y_velocity',                    &
                                                  long_name='y velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VEL) = netcdf_field_info(name='z_velocity',                    &
                                                  long_name='z velocity component',     &
                                                  std_name='',                          &
                                                  unit='m/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_X_VTEND) = netcdf_field_info(name='x_vorticity_tendency',        &
                                                    long_name='x vorticity tendency',   &
                                                    std_name='',                        &
                                                    unit='1/s',                         &
                                                    dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VTEND) = netcdf_field_info(name='y_vorticity_tendency',        &
                                                    long_name='y vorticity tendency',   &
                                                    std_name='',                        &
                                                    unit='1/s',                         &
                                                    dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VTEND) = netcdf_field_info(name='z_vorticity_tendency',        &
                                                    long_name='z vorticity tendency',   &
                                                    std_name='',                        &
                                                    unit='1/s',                         &
                                                    dtype=NF90_DOUBLE)

            nc_dset(NC_NPARG) = netcdf_field_info(name='nparg',                            &
                                                  long_name='number of parcels per cell',  &
                                                  std_name='',                             &
                                                  unit='1',                                &
                                                  dtype=NF90_INT)

            nc_dset(NC_NSPARG) = netcdf_field_info(name='nsparg',                                &
                                                   long_name='number of small parcels per cell', &
                                                   std_name='',                                  &
                                                   unit='1',                                     &
                                                   dtype=NF90_INT)

            nc_dset(NC_X_VOR) = netcdf_field_info(name='x_vorticity',                   &
                                                  long_name='x vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VOR) = netcdf_field_info(name='y_vorticity',                   &
                                                  long_name='y vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VOR) = netcdf_field_info(name='z_vorticity',                   &
                                                  long_name='z vorticity component',    &
                                                  std_name='',                          &
                                                  unit='1/s',                           &
                                                  dtype=NF90_DOUBLE)

            nc_dset(NC_TBUOY) = netcdf_field_info(name='buoyancy',                      &
                                                   long_name='total buoyancy',          &
                                                   std_name='',                         &
                                                   unit='m/s^2',                        &
                                                   dtype=NF90_DOUBLE)

#ifndef ENABLE_DRY_MODE
            nc_dset(NC_DBUOY) = netcdf_field_info(name='dry_buoyancy',                  &
                                                   long_name='dry buoyancy',            &
                                                   std_name='',                         &
                                                   unit='m/s^2',                        &
                                                   dtype=NF90_DOUBLE)

            nc_dset(NC_HUM) = netcdf_field_info(name='humidity',                        &
                                                long_name='specific humidity',          &
                                                std_name='',                            &
                                                unit='kg/kg',                           &
                                                dtype=NF90_DOUBLE   )

            nc_dset(NC_LBUOY) = netcdf_field_info(name='liquid_water_content',          &
                                                  long_name='liquid-water content',     &
                                                  std_name='',                          &
                                                  unit='1',                             &
                                                  dtype=NF90_DOUBLE)
#endif

            nc_dset(NC_VOL) = netcdf_field_info(name='volume',                          &
                                                long_name='volume',                     &
                                                std_name='',                            &
                                                unit='m^3',                             &
                                                dtype=NF90_DOUBLE)
        end subroutine set_netcdf_field_info

end module field_netcdf
