! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel_info
    use constants, only : zero, two
    use parcel_container, only : parcels, n_parcels
    use ellipse, only : get_ab
    use parameters, only : dx, vcell, ncell, extent, lower, nx, nz
    implicit none

    private :: init_random_positions,  &
               init_regular_positions

    contains


        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine parcel_default
            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel_info%n_per_cell * ncell

            if (parcel_info%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif

            ! initialize the volume of each parcel
            parcels%volume(1:n_parcels, 1) = vcell / dble(parcel_info%n_per_cell)

            if (parcel_info%is_elliptic) then
                deallocate(parcels%stretch)


                ! initialze circles
                parcels%B(1:n_parcels, 1) = get_ab(parcels%volume(1:n_parcels, 1)) ! B11
                parcels%B(1:n_parcels, 2) = zero                                   ! B12

            else
                deallocate(parcels%B)
                parcels%stretch = zero
            endif

            parcels%velocity(1:n_parcels, :) = zero
            parcels%buoyancy(1:n_parcels, 1) = zero
            parcels%humidity(1:n_parcels, 1) = zero
        end subroutine parcel_default


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
            n_per_dim = dsqrt(dble(parcel_info%n_per_cell))
            if (n_per_dim ** 2 .ne. parcel_info%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel_info%n_per_cell, ") not a square."
                stop
            endif

            del = dx / dble(two * n_per_dim)

            k = 1
            do j = 0, nz-1
                do i = 0, nx-1
                    do jj = 1, 2 * n_per_dim, 2
                        do ii = 1, 2 * n_per_dim, 2
                            parcels%position(k, 1) = lower(1) + dble(i) * dx(1) + del(1) * dble(ii)
                            parcels%position(k, 2) = lower(2) + dble(j) * dx(2) + del(2) * dble(jj)
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo
        end subroutine init_regular_positions

end module parcel_init
