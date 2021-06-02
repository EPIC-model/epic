! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel_info
    use constants, only : zero, two, one
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parameters, only : dx, vcell, ncell, extent, lower, nx, nz
    implicit none

    private :: init_random_positions,  &
               init_regular_positions

    contains


        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine parcel_default
            double precision :: lam, ratio

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel_info%n_per_cell * ncell

            if (parcel_info%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif

            ! initialize the volume of each parcel
            parcels%volume(1:n_parcels) = vcell / dble(parcel_info%n_per_cell)

            if (parcel_info%is_elliptic) then
                deallocate(parcels%stretch)

                ratio = dx(1) / dx(2)

                ! aspect ratio: lam = a / b
                lam = max(dx(2) / dx(1), ratio)

                ! B11
                parcels%B(1:n_parcels, 1) = ratio * get_ab(parcels%volume(1:n_parcels))

                ! B12
                parcels%B(1:n_parcels, 2) = zero

                call init_refine(lam)

            else
                deallocate(parcels%B)
                parcels%stretch = zero
            endif

            parcels%velocity(1:n_parcels, :) = zero
            parcels%vorticity(1:n_parcels, :) = zero
            parcels%buoyancy(1:n_parcels) = zero
            parcels%humidity(1:n_parcels) = zero
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

        subroutine init_refine(lam)
            double precision, intent(inout) :: lam
            double precision                :: B22, a2

            if (.not. parcel_info%is_elliptic) then
                return
            endif

            ! do refining by splitting
            do while (lam >= parcel_info%lambda)
                call split_ellipses(parcels, parcel_info%lambda, parcel_info%vmaxfraction)
                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine

end module parcel_init
