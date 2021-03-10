module parser
    use types, only : mesh_type, time_info_type, parcel_info_type
    use parameters, only : filename, h5freq
    implicit none

    type(mesh_type) :: mesh

    ! declare some subroutine private
    private :: parse_domain_info,  &
               parse_output_info,  &
               parse_parcel_info,  &
               parse_stepper_info

    contains
        ! Parse OUTPUT_INFO namelist
        subroutine parse_output_info(ios, fn)
            integer, intent(inout)   :: ios
            integer, intent(in)      :: fn

            ! Namelist definition
            namelist /OUTPUT_INFO/ h5freq

            read(nml=OUTPUT_INFO, iostat=ios, unit=fn)
        end subroutine parse_output_info

        ! Parse DOMAIN_INFO namelist
        subroutine parse_domain_info(ios, fn)
            integer,            intent(inout)   :: ios
            integer,            intent(in)      :: fn
            double precision                    :: origin(2)
            double precision                    :: extent(2)
            integer                             :: grid(2)
            character(len=16)                   :: bc(2)

            namelist /DOMAIN_INFO/ origin, extent, grid, bc

            read(nml=DOMAIN_INFO, iostat=ios, unit=fn)

            mesh%origin = origin
            mesh%extent = extent
            mesh%grid = grid
            mesh%bc = bc
        end subroutine

        ! Parse PARCEL_INFO namelist
        subroutine parse_parcel_info(ios, fn, parcel)
            integer,                intent(inout)   :: ios
            integer,                intent(in)      :: fn
            type(parcel_info_type), intent(inout)   :: parcel

            namelist /PARCEL_INFO/ parcel

            read(nml=PARCEL_INFO, iostat=ios, unit=fn)
        end subroutine

        ! Parse STEPPER_INFO namelist
        subroutine parse_stepper_info(ios, fn, time)
            integer,              intent(inout) :: ios
            integer,              intent(in)    :: fn
            type(time_info_type), intent(inout) :: time

            namelist /STEPPER_INFO/ time

            read(nml=STEPPER_INFO, iostat=ios, unit=fn)
        end subroutine


        ! Parsing namelists
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file(time, parcel)
            integer                                :: ios
            integer                                :: fn = 1
            type(time_info_type),   intent(inout)  :: time
            type(parcel_info_type), intent(inout)  :: parcel

            ! Check whether file exists.
            inquire(file=filename, iostat=ios)

            if (ios /= 0) then
                write (*, *) 'Error: input file "', trim(filename), '" does not exist.'
                return
            endif

            ! Open and read Namelist file.
            open (action='read', file=filename, iostat=ios, newunit=fn)

            call parse_output_info(ios, fn)
            call parse_domain_info(ios, fn)
            call parse_parcel_info(ios, fn, parcel)
            call parse_stepper_info(ios, fn, time)

            if (ios /= 0) then
                write (*, *) 'Error: invalid Namelist format.'
            end if

            close(fn)

        end subroutine read_config_file


end module parser
