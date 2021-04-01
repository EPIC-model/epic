module parcel_bc
    use constants, only : lower, upper, extent
    use parcel_container, only : n_parcels
    implicit none

    contains

        ! apply boundary condition in all directions
        ! x -- periodic bc
        ! y -- free slip bc
        subroutine apply_parcel_bc(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            integer                         :: n

            do n = 1, n_parcels
                ! horizontal bc
                call do_periodic(position(n, :))

                ! vertical bc
                call do_free_slip(position(n, :), velocity(n, :))
            enddo
        end subroutine apply_parcel_bc


        ! apply periodic boundary condition in x (horizontal) direction
        subroutine apply_periodic_bc(position)
            double precision, intent(inout) :: position(:, :)
            integer                         :: n

            do n = 1, n_parcels
                call do_periodic(position(n, :))
            enddo
        end subroutine apply_periodic_bc


        ! apply periodic boundary condition in y (vertical) direction
        subroutine apply_free_slip_bc(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            integer                         :: n

            do n = 1, n_parcels
                call do_free_slip(position(n, :), velocity(n, :))
            enddo

        end subroutine apply_free_slip_bc

        ! apply periodic bc on n-th parcel
        subroutine do_periodic(position)
            double precision, intent(inout) :: position(2)

            if (position(1) > upper(1)) then
                position(1) = position(1) - extent(1)
            endif

            if (position(1) < lower(1)) then
                position(1) = position(1) + extent(1)
            endif
        end subroutine do_periodic


        ! apply free slip bc on n-th parcel
        subroutine do_free_slip(position, velocity)
            double precision, intent(inout) :: position(2), velocity(2)

            if (position(2) >= upper(2)) then
                velocity(2) = 0.0
                position(2) = 2.0 * upper(2) - position(2)
            endif

            if (position(2) <= lower(2)) then
                velocity(2) = 0.0
                position(2) = 2.0 * lower(2) - position(2)
            endif
        end subroutine do_free_slip

end module parcel_bc
