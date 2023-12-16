! =============================================================================
!                       Parcel boundary conditions
!                       periodic in x (horizontal)
!                       reflective in z (vertical)
! =============================================================================
module parcel_bc
    use constants, only : zero, two
    use parameters, only : lower, upper, extent, hli, center
    use parcel_container, only : n_parcels, parcels
    use surface_parcel_bc, only : apply_surface_parcel_bc
    use omp_lib
    implicit none

    contains

        ! Apply periodic bc on n-th parcel (horizontally)
        ! @param[inout] position vector of parcel
        subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(2)
            position(1) = position(1) - extent(1) * dble(int((position(1) - center(1)) * hli(1)))
        end subroutine apply_periodic_bc

        ! Apply mirroring bc on n-th parcel (vertically)
        ! @param[inout] position vector of parcel
        ! @param[inout] B matrix of parcel
        subroutine apply_reflective_bc(position, B)
            double precision, intent(inout) :: position(2), B(2)

            if (position(2) > upper(2)) then
                position(2) = two * upper(2) - position(2)
                B(2) = -B(2)
            else if (position(2) < lower(2)) then
                position(2) = two * lower(2) - position(2)
                B(2) = -B(2)
            endif
        end subroutine apply_reflective_bc

        ! Apply all boundary conditions to all parcels
        subroutine apply_parcel_bc
            integer :: n

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! horizontal direction
                call apply_periodic_bc(parcels%position(:, n))

                ! vertical direction
                call apply_reflective_bc(parcels%position(:, n), parcels%B(:, n))
            enddo
            !$omp end do
            !$omp end parallel

            call apply_surface_parcel_bc

        end subroutine apply_parcel_bc

end module parcel_bc
