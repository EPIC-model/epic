program epic
    use model, only : pre_run, run, post_run
    use parameters, only : filename, verbose
    implicit none

    ! Read command line (verbose, filename, etc.)
    call parse_command_line(filename, verbose)

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

end program epic


! Get the file name provided via the command line
subroutine parse_command_line(filename, verbose)
    integer                          :: i
    character(len=32)                :: arg
    character(len=32), intent(inout) :: filename
    logical,           intent(inout) :: verbose

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

