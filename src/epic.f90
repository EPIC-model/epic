program epic
    use model, only : pre_run, run, post_run
    implicit none

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

end program epic


! Get the file name provided via the command line
subroutine parse_command_line
    use parameters, only : filename, verbose
    integer                          :: i
    character(len=32)                :: arg

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
        print *, 'No configuration file provided. Run code with "./epic --config [config file]"'
        stop
    endif

    ! This is the main application of EPIC
    if (verbose) then
        print *, 'Running EPIC with "', filename, '"'
    endif
end subroutine parse_command_line

