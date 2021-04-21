! =============================================================================
!               This module initializes all parcels and fields.
! =============================================================================
module init
    use options, only : parcel_info, time
    use parameters, only : dx, vcell, ncell, extent, lower, nx, nz
    use fields, only : velocity_f,      &
                       strain_f,        &
                       volg,            &
                       vorticity_f,     &
                       get_position
    use ellipse, only : get_ab
    use parcel_container, only : parcels, n_parcels
    implicit none

    private :: init_random_positions,  &
               init_regular_positions, &
               init_stretch,           &
               init_B_matrix,          &
               init_velocity,          &
               init_velocity_field

    contains

        !
        ! parcel initialize subroutines
        !

        subroutine init_parcels
            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel_info%n_per_cell * ncell

            if (parcel_info%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif

            call init_stretch

            ! initialize the volume of each parcel
            parcels%volume(1:n_parcels, 1) = vcell / parcel_info%n_per_cell

            call init_B_matrix

            call init_velocity


        end subroutine init_parcels


        subroutine init_random_positions
            double precision :: val
            integer          :: n, k
            integer, allocatable :: seed(:)

            call random_seed(size=k)
            allocate(seed(1:k))
            seed(:) = parcel_info%seed
            call random_seed(put=seed)

            do n = 1, n_parcels
                call random_number(val)
                parcels%position(n, 1)= lower(1) + val * extent(1)
                call random_number(val)
                parcels%position(n, 2) = lower(2) + val * extent(2)
            enddo

            deallocate(seed)
        end subroutine init_random_positions

        subroutine init_regular_positions
            integer          :: i, ii, j, jj, k, n_per_dim
            double precision :: del(2)

            ! number of parcels per dimension
            n_per_dim = sqrt(parcel_info%n_per_cell + 0.0)
            if (n_per_dim ** 2 .ne. parcel_info%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel_info%n_per_cell, ") not a square."
                stop
            endif

            del = dx / (2.0 * n_per_dim)

            k = 1
            do j = 0, nz-1
                do i = 0, nx-1
                    do jj = 1, 2 * n_per_dim, 2
                        do ii = 1, 2 * n_per_dim, 2
                            parcels%position(k, 1) = lower(1) + i * dx(1) + del(1) * ii
                            parcels%position(k, 2) = lower(2) + j * dx(2) + del(2) * jj
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
                ! initialze circles
                parcels%B(1:n_parcels, 1) = get_ab(parcels%volume(1:n_parcels, 1)) ! B11
                parcels%B(1:n_parcels, 2) = 0.0                                    ! B12
            else
                deallocate(parcels%B)
            endif
        end subroutine init_B_matrix

        subroutine init_velocity
            use taylorgreen, only : get_flow_velocity
            integer :: n

            do n = 1, n_parcels
                parcels%velocity(n, :) = get_flow_velocity(parcels%position(n, :))
            enddo
        end subroutine init_velocity

        !
        ! field initialize subroutines
        !

        subroutine init_fields
            call init_velocity_field

            allocate(volg(-1:nz+1, 0:nx-1, 1))
            volg = 0.0

            if (time%is_adaptive) then
                call init_vorticity_field
            endif

        end subroutine init_fields

        subroutine init_velocity_field
            use taylorgreen, only : get_flow_velocity, &
                                    get_flow_gradient
            integer :: i, j
            double precision :: pos(2)

            allocate(velocity_f(0:nx-1, -1:nz+1, 2))

            allocate(strain_f(0:nx-1, -1:nz+1, 4))

            do j = -1, nz+1
                do i = 0, nx-1
                    pos = get_position((/i, j/))

                    velocity_f(i, j, :) = get_flow_velocity(pos)

                    strain_f(i, j, :) = get_flow_gradient(pos)
                enddo
            enddo
        end subroutine init_velocity_field

        subroutine init_vorticity_field
            use taylorgreen, only : get_flow_vorticity
            integer :: i, j
            double precision :: pos(2)

            ! vorticity has no halo grid points in y
            allocate(vorticity_f(0:nx-1, 0:nz, 1))

            do j = 0, nz
                do i = 0, nx-1
                    pos = get_position((/i, j/))

                    vorticity_f(i, j, :) = get_flow_vorticity(pos)
                enddo
            enddo

        end subroutine init_vorticity_field

end module init
