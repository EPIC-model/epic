! =============================================================================
!               This program writes fields to HDF5 in EPIC format.
! =============================================================================
program models
    use taylorgreen
    use straka
    use robert
    implicit none

    character(len=32) :: filename = ''
    logical           :: verbose = .false.

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    call generate_fields(trim(filename))

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
