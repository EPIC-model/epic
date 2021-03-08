module parser
    implicit none

    logical :: verbose = .false.
    character(:), allocatable :: filename

    integer :: h5freq

    type :: coord_type
        double precision :: x
        double precision :: y
    end type coord_type

    type :: aint_type
        integer :: x
        integer :: y
    end type aint_type

    type :: astring_type
        character(len=32) :: x
        character(len=32) :: y
    end type astring_type

    type(coord_type) :: origin

    type(coord_type) :: extent

    type(aint_type) :: grid

    type(astring_type) :: bc



    contains
        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer :: i
            character(len=32) :: arg

            filename = ''
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
                endif
                i = i+1
            end do

            if (filename == '') then
                write(*,*) 'No configuration file provided. Run code with "./epic --config [config file]"'
                return
            endif

            ! This is the main application of EPIC
            if (verbose) then
                print *, 'Running EPIC with "', filename, '"'
            endif
        end subroutine parse_command_line

        ! Parse output namelist
        subroutine get_output(ios, fn)
            integer, intent(inout)   :: ios
            integer, intent(in)      :: fn

            ! Namelist definition
            namelist /output/ h5freq

            read(nml=output, iostat=ios, unit=fn)
        end subroutine get_output

        ! Parse domain namelist
        subroutine get_domain(ios, fn)
            integer,            intent(inout)   :: ios
            integer,            intent(in)      :: fn

            namelist /domain/ origin, extent, grid, bc

            read(nml=domain, iostat=ios, unit=fn)
        end subroutine


        ! Parsing namelists
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1

            ! set all parser attributes (verbose, filename, etc.)
            call parse_command_line()


            ! Check whether file exists.
            inquire(file=filename, iostat=ios)

            if (ios /= 0) then
                write (*, *) 'Error: input file "', trim(filename), '" does not exist.'
                return
            endif

            ! Open and read Namelist file.
            open (action='read', file=filename, iostat=ios, newunit=fn)

!             call get_output(ios, fn)
            call get_domain(ios, fn)

            if (ios /= 0) then
                write (*, *) 'Error: invalid Namelist format.'
            end if

            close(fn)

        end subroutine read_config_file


end module parser
