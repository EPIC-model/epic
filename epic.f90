program epic
    use parser, only : get_config_file
    implicit none

    ! Obtain the file provided via the command line
    character(:), allocatable :: filename
    call get_config_file(filename)

    ! This is the main application of EPIC
    print *, 'Running EPIC with "', filename, '"'
end program epic
