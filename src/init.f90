module init
    use parcels
    implicit none

    private :: init_random_positions, &
               init_regular_positions

    contains

        subroutine init_parcels
            use options, only : is_random

            if (is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif
        end subroutine init_parcels


        subroutine init_random_positions

        end subroutine init_random_positions


        subroutine init_regular_positions

        end subroutine init_regular_positions

end module init
