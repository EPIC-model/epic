! =============================================================================
!                       Parcel boundary conditions
!                       periodic in x (horizontal)
!                       free slip in z (vertical)
! =============================================================================
module parcel_bc
    use constants, only : zero, two
    use parameters, only : lower, upper, extent, hli, center
    use parcel_container, only : n_parcels
    use omp_lib
    implicit none

    contains

        ! apply boundary condition in all directions
        ! x -- periodic bc
        ! z -- free slip bc
        subroutine apply_parcel_bc(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            integer                         :: n

            !$omp parallel num_threads(4)
            !$omp do private(n)
            do n = 1, n_parcels
                ! horizontal bc
                call apply_periodic_bc(position(n, :))

                ! vertical bc
                call apply_free_slip_bc(position(n, :), velocity(n, :))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine apply_parcel_bc


        ! apply periodic bc on n-th parcel
        subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(2)
            position(1) = position(1) - extent(1) * dble(int((position(1) - center(1)) * hli(1)))
        end subroutine apply_periodic_bc


        ! apply free slip bc on n-th parcel
        subroutine apply_free_slip_bc(position, velocity)
            double precision, intent(inout) :: position(2), velocity(2)

            if (position(2) >= upper(2)) then
                velocity(2) = zero
                position(2) = two * upper(2) - position(2)
            endif

            if (position(2) <= lower(2)) then
                velocity(2) = zero
                position(2) = two * lower(2) - position(2)
            endif
        end subroutine apply_free_slip_bc


        subroutine apply_reflective_bc(position, B)
            double precision, intent(inout) :: position(2), B(2)

            if (position(2) > upper(2)) then
                position(2) = two * upper(2) - position(2)
                B(2) = -B(2)
            endif

            if (position(2) < lower(2)) then
                position(2) = two * lower(2) - position(2)
                B(2) = -B(2)
            endif
        end subroutine apply_reflective_bc

end module parcel_bc
