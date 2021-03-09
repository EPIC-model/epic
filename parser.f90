module parser
    use types, only : mesh_type
    implicit none

    logical :: verbose = .false.
    character(:), allocatable :: filename

    integer :: h5freq

    type(mesh_type) :: mesh

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
            double precision                    :: origin(2)
            double precision                    :: extent(2)
            integer                             :: grid(2)
            character(len=16)                   :: bc(2)

            namelist /domain/ origin, extent, grid, bc

            read(nml=domain, iostat=ios, unit=fn)

            mesh%origin = origin
            mesh%extent = extent
            mesh%grid = grid
            mesh%bc = bc
        end subroutine

        ! Parse parcel namelist
        subroutine get_parcel(ios, fn)
            integer,            intent(inout)   :: ios
            integer,            intent(in)      :: fn
            integer                             :: n_per_cell = 4
            logical                             :: random     = .false.
            integer                             :: seed       = 42
            logical                             :: elliptic   = .true.

            namelist /parcel/ n_per_cell, random, seed, elliptic

            read(nml=parcel, iostat=ios, unit=fn)

            write(*,*) elliptic

        end subroutine

        ! Parse time namelist
        subroutine get_time(ios, fn)
            integer,            intent(inout)   :: ios
            integer,            intent(in)      :: fn
            double precision                    :: limit    = 0.0
            double precision                    :: max_step = 0.0
            logical                             :: adaptive = .true.

            namelist /time/ limit, max_step, adaptive

            read(nml=time, iostat=ios, unit=fn)
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

            call get_output(ios, fn)
            call get_domain(ios, fn)
            call get_parcel(ios, fn)
            call get_time(ios, fn)

            if (ios /= 0) then
                write (*, *) 'Error: invalid Namelist format.'
            end if

            close(fn)

        end subroutine read_config_file


end module parser
