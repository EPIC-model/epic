! =============================================================================
!               This program writes fields to HDF5 in EPIC format.
! =============================================================================
program models
    use taylorgreen
    use straka
    use robert
    implicit none

    character(len=32) :: filename = ''
    character(len=32) :: model = ''
    character(len=32) :: h5fname = ''
    logical           :: verbose = .false.

    type box_type
        integer          :: nc(2)       ! number of cells
        double precision :: extent(2)   ! size of domain
        double precision :: origin(2)   ! origin of domain (lower left corner)
    end type box_type

    type(box_type) :: box


    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    call read_config_file

    call generate_fields(trim(model))

    contains

        subroutine generate_fields(name)
            character(*), intent(in) :: name

            select case (trim(name))
                case ('TaylorGreen')
                    call taylorgreen_init
                case ('Straka')
                    call straka_init
                case ('Robert')
                    call robert_init
                case default
                    print *, "Invalid simulation type: '", trim(name), "'"
                    stop
            end select
        end subroutine generate_fields


        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /MODELS/ model, h5fname, box

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

            ! check whether h5 file already exists
            inquire(file=h5fname, exist=exists)

            if (exists) then
                print *, 'Error: output file "', trim(h5fname), '" already exists.'
                stop
            endif
        end subroutine read_config_file

        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer                          :: i
            character(len=32)                :: arg

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
                print *, 'No configuration file provided. Run code with "./models --config [config file]"'
                stop
            endif
        end subroutine parse_command_line
end program models
