module init
    use parameters, only : mesh
    use parcels
    implicit none

    private :: init_random_positions, &
               init_regular_positions

    contains

        subroutine init_parcels
            use parameters, only : is_random
            integer :: n_cells

            n_cells = product(mesh%grid)

            ! set the number of parcels (see parcels.f90)
            ! we use 4 parcels per grid cell
            n_parcels = 4 * n_cells

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

            do i = 1, n_parcels
                call random_number(val)
                x(i)= mesh%origin(1) + val
                call random_number(val)
                y(i) = mesh%origin(2) + val
            enddo
        end subroutine init_random_positions

        subroutine init_regular_positions
            integer :: i, j, k, l
            double precision :: dx(2)

            dx = mesh%extent / (mesh%grid - 1)

            k = 1
            do i = 1, mesh%grid(1)
                do j = 1, mesh%grid(2)
                    l = mod(k, 2)
                    x(k) = mesh%origin(1) + (0.25 + i + 0.5 * l) * dx(1)
                    y(k) = mesh%origin(2) + (0.25 + j + 0.5 * l) * dx(2)
                    k = k + 1
                enddo
            enddo

        end subroutine init_regular_positions

end module init
