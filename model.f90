module model
    use parser, only : get_config_file
    implicit none

    contains
        subroutine setup
            ! Obtain the file provided via the command line
            character(:), allocatable :: filename
            call get_config_file(filename)


        end subroutine


        subroutine run

        end subroutine run


end module model
