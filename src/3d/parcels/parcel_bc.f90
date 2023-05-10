! =============================================================================
!                       Parcel boundary conditions
!                       periodic in x (zonal) and in y (meridional)
!                       reflective in z (vertical)
! =============================================================================
module parcel_bc
    use constants, only : zero, two
    use parameters, only : lower, upper, extent, hli, center
    use parcel_container, only : n_parcels, parcels
    use parcel_mpi, only : parcel_communicate
    use omp_lib
    implicit none

    contains

        ! Apply periodic bc on n-th parcel (zonally and meridionally)
        ! @param[inout] position vector of parcel
        pure subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(3)
            position(1) = position(1) - extent(1) * dble(int((position(1) - center(1)) * hli(1)))
            position(2) = position(2) - extent(2) * dble(int((position(2) - center(2)) * hli(2)))
        end subroutine apply_periodic_bc

        ! Apply mirroring bc on n-th parcel (vertically)
        ! @param[inout] position vector of parcel
        ! @param[inout] B matrix of parcel
        pure subroutine apply_reflective_bc(position, B)
            double precision, intent(inout) :: position(3), B(5)

            if (position(3) > upper(3)) then
                position(3) = two * upper(3) - position(3)
                ! flip sign of B13 and B23
                B(3) = -B(3)
                B(5) = -B(5)
            else if (position(3) < lower(3)) then
                position(3) = two * lower(3) - position(3)
                ! flip sign of B13 and B23
                B(3) = -B(3)
                B(5) = -B(5)
            endif
        end subroutine apply_reflective_bc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Apply all boundary conditions to all parcels
        subroutine apply_parcel_bc
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! zonal direction
                call apply_periodic_bc(parcels%position(:, n))

                ! vertical direction
                call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_parcel_bc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Apply periodic boundary condition to all parcels
        ! @param[inout] position vector of parcels
        subroutine apply_parcel_periodic_bc
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! zonal and meridional direction
                call apply_periodic_bc(parcels%position(:, n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_parcel_periodic_bc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Performs a parcel swap and then applies the periodic boundary condition
        ! to all parcels
        ! @param[in] pindex are the parcel indices to check (optional)
        subroutine apply_swap_periodicity(pindex)
            integer, optional, intent(in) :: pindex(:)
            integer                       :: n

            call parcel_communicate(pindex)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! zonal and meridional direction
                call apply_periodic_bc(parcels%position(:, n))
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_swap_periodicity

end module parcel_bc
