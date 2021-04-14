! =============================================================================
!                       Parcel boundary conditions
!                       periodic in x (horizontal)
!                       free slip in z (vertical)
! =============================================================================
module parcel_bc
    use parameters, only : lower, upper, extent, hli
    use parcel_container, only : n_parcels, are_parcels_modified
    implicit none

    contains

        ! apply boundary condition in all directions
        ! x -- periodic bc
        ! z -- free slip bc
        subroutine apply_parcel_bc(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            integer                         :: n

            do n = 1, n_parcels
                ! horizontal bc
                call apply_periodic_bc(position(n, :))

                ! vertical bc
                call apply_free_slip_bc(position(n, :), velocity(n, :))
            enddo

            are_parcels_modified = .true.

        end subroutine apply_parcel_bc


        ! apply periodic bc on n-th parcel
        subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(2)
            position(1) = position(1) - extent(1) * dble(int(position(1) * hli(1)))
        end subroutine apply_periodic_bc


        ! apply free slip bc on n-th parcel
        subroutine apply_free_slip_bc(position, velocity)
            double precision, intent(inout) :: position(2), velocity(2)

            if (position(2) >= upper(2)) then
                velocity(2) = 0.0
                position(2) = 2.0 * upper(2) - position(2)
            endif

            if (position(2) <= lower(2)) then
                velocity(2) = 0.0
                position(2) = 2.0 * lower(2) - position(2)
            endif
        end subroutine apply_free_slip_bc

end module parcel_bc
