module init
    use parameters, only : mesh
    use parcels
    implicit none

    private :: init_random_positions, &
               init_regular_positions

    contains

        subroutine init_parcels
            use parameters, only : is_random, n_per_cell
            integer :: n_cells

            n_cells = product(mesh%grid)

            print *, n_cells

            n_parcels = n_cells * n_per_cell

            if (is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif
        end subroutine init_parcels


        subroutine init_random_positions
            use parameters, only : seed
            double precision :: val
            integer :: i

            call random_seed !put=seed)

            do i=1, n_parcels
                call random_number(val)
                x(i)= mesh%origin(1) + val
                call random_number(val)
                y(i) = mesh%origin(2) + val
            enddo
        end subroutine init_random_positions


        subroutine init_regular_positions

        end subroutine init_regular_positions

end module init
