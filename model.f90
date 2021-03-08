module model
    use parser, only : read_config_file
    implicit none

    contains
        subroutine setup
            ! parse the config file
            call read_config_file()


        end subroutine


        subroutine run

        end subroutine run


end module model
