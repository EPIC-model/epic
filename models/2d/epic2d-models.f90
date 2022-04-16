! =============================================================================
!               This program writes fields to NetCDF in EPIC format.
! =============================================================================
program epic2d_models
    use taylor_green_2d
    use straka_2d
    use robert_2d
    use constants, only : pi, zero
    use parameters, only : nx, nz, dx, lower, extent
    use netcdf_utils
    use netcdf_writer
    use config, only : package_version, cf_version
    use physics, only : read_physical_quantities_from_namelist
    implicit none

    logical            :: verbose = .false.
    character(len=512) :: filename = ''
    character(len=512) :: model = ''
    character(len=512) :: ncfname = ''
    integer            :: ncid
    integer            :: dimids(3), axids(3)

    type box_type
        integer          :: ncells(2)   ! number of cells
        double precision :: extent(2)   ! size of domain
        double precision :: origin(2)   ! origin of domain (lower left corner)
    end type box_type

    type(box_type) :: box


    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    call read_config_file

    call read_physical_quantities_from_namelist(trim(filename))

    call generate_fields

    contains

        subroutine generate_fields

            call create_netcdf_file(ncfname, .false., NF90_NETCDF4, ncid)

            dx = box%extent / dble(box%ncells)
            nx = box%ncells(1)
            nz = box%ncells(2)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   epic_version=package_version, &
                                   file_type='fields',           &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, box%ncells)

            call define_netcdf_spatial_dimensions_2d(ncid=ncid,            &
                                                     ncells=box%ncells,    &
                                                     dimids=dimids(1:2),   &
                                                     axids=axids(1:2))

            call define_netcdf_temporal_dimension(ncid, dimids(3), axids(3))

            if (model == 'TaylorGreen') then
                ! make origin and extent always a multiple of pi
                box%origin = pi * box%origin
                box%extent = pi * box%extent
                dx = dx * pi
            endif

            ! write box
            lower = box%origin
            extent = box%extent
            call write_netcdf_box(ncid, lower, extent, box%ncells)

            select case (trim(model))
                case ('TaylorGreen')
                    call taylor_green_init(ncid, dimids, nx, nz, box%origin, dx)
                case ('Straka')
                    call straka_init(ncid, dimids, nx, nz, box%origin, dx)
                case ('Robert')
                    call robert_init(ncid, dimids, nx, nz, box%origin, dx)
                case default
                    print *, "Unknown model: '", trim(model), "'."
                    stop
            end select

            call write_netcdf_axis_2d(ncid, dimids, lower, dx, box%ncells)

            ! write time
            call write_netcdf_scalar(ncid, axids(3), zero, 1)

            call close_netcdf_file(ncid)
        end subroutine generate_fields


        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /MODELS/ model, ncfname, box, tg_flow, straka_flow, robert_flow

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                print *, 'Error: input file "', trim(filename), '" does not exist.'
                stop
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=MODELS, iostat=ios, unit=fn)

            if (ios /= 0) then
                print *, 'Error: invalid Namelist format.'
                stop
            end if

            close(fn)

            ! check whether NetCDF file already exists
            inquire(file=ncfname, exist=exists)

            if (exists) then
                print *, 'Error: output file "', trim(ncfname), '" already exists.'
                stop
            endif
        end subroutine read_config_file

        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer                          :: i
            character(len=512)               :: arg

            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--config') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    filename = trim(arg)
                else if (arg == '--verbose') then
                    verbose = .true.
                else if (arg == '--help') then
                    print *, 'Run code with "epic2d-models --config [config file]"'
                    stop
                endif
                i = i+1
            end do

            if (filename == '') then
                print *, 'No configuration file provided. Run code with "epic2d-models --config [config file]"'
                stop
            endif
        end subroutine parse_command_line
end program epic2d_models
