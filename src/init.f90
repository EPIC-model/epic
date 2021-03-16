module init
    use parameters, only : mesh, parcel_info
    use fields, only : velocity, get_mesh_spacing
    use parcel_container, only : parcels, n_parcels
    implicit none

    private :: init_random_positions,  &
               init_regular_positions, &
               init_stretch,           &
               init_B_matrix,          &
               init_velocity

    contains

        !
        ! parcel initialize subroutines
        !

        subroutine init_parcels
            integer          :: n_cells
            double precision :: cell_volume

            n_cells = product(mesh%grid - 1)

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel_info%n_per_cell * n_cells

            if (parcel_info%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif

            call init_stretch

            call init_B_matrix

            call init_velocity

            ! initialize the volume of each parcel
            cell_volume = product(get_mesh_spacing())
            parcels%volume = cell_volume / parcel_info%n_per_cell

        end subroutine init_parcels


        subroutine init_random_positions
            double precision :: val
            integer :: i

            call random_seed !put=parcel_info%seed)

            do i = 1, n_parcels
                call random_number(val)
                parcels%position(i, 1)= mesh%origin(1) + val
                call random_number(val)
                parcels%position(i, 2) = mesh%origin(2) + val
            enddo
        end subroutine init_random_positions

        subroutine init_regular_positions
            integer          :: i, ii, j, jj, k, n_per_dim
            double precision :: dx(2), del(2)

            dx = get_mesh_spacing()

            ! number of parcels per dimension
            n_per_dim = sqrt(parcel_info%n_per_cell + 0.0)
            if (n_per_dim ** 2 .ne. parcel_info%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel_info%n_per_cell, ") not a square."
                stop
            endif

            del = dx / (n_per_dim + 1)

            k = 1
            do j = 0, mesh%grid(2) - 2
                do i = 0, mesh%grid(1) - 2
                    do jj = 1, n_per_dim
                        do ii = 1, n_per_dim
                            parcels%position(k, 1) = mesh%origin(1) + i * dx(1) + del(1) * ii
                            parcels%position(k, 2) = mesh%origin(2) + j * dx(2) + del(2) * jj
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo
        end subroutine init_regular_positions

        subroutine init_stretch
            if (parcel_info%is_elliptic) then
                deallocate(parcels%stretch)
            else
                parcels%stretch = 0.0
            endif
        end subroutine init_stretch

        subroutine init_B_matrix
            if (parcel_info%is_elliptic) then
                parcels%B11 = 1.0
                parcels%B12 = 0.0
            else
                deallocate(parcels%B11)
                deallocate(parcels%B12)
            endif
        end subroutine init_B_matrix

        subroutine init_velocity
            use taylorgreen, only : get_flow_velocity
            integer :: i
            double precision :: dx(2), v(2)

            dx = get_mesh_spacing()

            do i = 1, n_parcels
                parcels%velocity(i, :) = get_flow_velocity(parcels%position(i, :))
            enddo
        end subroutine init_velocity

        !
        ! field initialize subroutines
        !

        subroutine init_velocity_field
            use taylorgreen, only : get_flow_velocity
            integer :: i, j
            double precision :: pos(2)
            double precision :: dx(2)

            allocate(velocity(mesh%grid(1), mesh%grid(2), 2))

            dx = get_mesh_spacing()

            do j = 1, mesh%grid(2)
                do i = 1, mesh%grid(1)
                    pos(1) = mesh%origin(1) + i * dx(1)
                    pos(2) = mesh%origin(2) + j * dx(2)
                    velocity(i, j, :) = get_flow_velocity(pos)
                enddo
            enddo
        end subroutine init_velocity_field

end module init
