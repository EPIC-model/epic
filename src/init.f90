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
            integer :: n_cells

            n_cells = product(mesh%grid - 1)

            ! set the number of parcels (see parcels.f90)
            ! we use 4 parcels per grid cell
            n_parcels = 4 * n_cells

            if (parcel_info%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif


            call init_stretch

            call init_velocity

        end subroutine init_parcels


        subroutine init_random_positions
            double precision :: val
            integer :: i

            call random_seed !put=parcel_info%seed)

            do i = 1, n_parcels
                call random_number(val)
                parcels%pos(i, 1)= mesh%origin(1) + val
                call random_number(val)
                parcels%pos(i, 2) = mesh%origin(2) + val
            enddo
        end subroutine init_random_positions

        subroutine init_regular_positions
            integer :: i, ii, j, jj k
            double precision :: dx(2), del(2)

            dx = get_mesh_spacing()

            del = dx / 3.0

            k = 1
            do j = 1, 2 * (mesh%grid(2) - 1)
                do i = 1, 2 * (mesh%grid(1) - 1)
                    ! cell shifts
                    ii = dx(1) * floor((i - 1) / 2)
                    jj = dx(2) * floor((j - 1) / 2)

                    parcels%pos(k, 1) = mesh%origin(1) + i * del(1) + ii
                    parcels%pos(k, 2) = mesh%origin(2) + j * del(2) + jj
                    k = k + 1
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
            double precision :: x, y
            double precision :: dx(2), v(2)

            dx = get_mesh_spacing()

            do i = 1, n_parcels
                parcels%vel(i, :) = get_flow_velocity(parcels%pos(i, 1), parcels%pos(i, 1))
            enddo
        end subroutine init_velocity

        !
        ! field initialize subroutines
        !

        subroutine init_velocity_field
            use taylorgreen, only : get_flow_velocity
            integer :: i, j
            double precision :: x, y
            double precision :: dx(2)

            allocate(velocity(mesh%grid(1), mesh%grid(2), 2))

            dx = get_mesh_spacing()

            do j = 1, mesh%grid(2)
                do i = 1, mesh%grid(1)
                    x = mesh%origin(1) + i * dx(1)
                    y = mesh%origin(2) + j * dx(2)
                    velocity(i, j, :) = get_flow_velocity(x, y)
                enddo
            enddo
        end subroutine init_velocity_field

end module init
