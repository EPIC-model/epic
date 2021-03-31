module parcel_bc
    use parameters, only : mesh
    use parcel_container, only : n_parcels
    implicit none

    private :: do_free_slip, &
               do_periodic

    contains

        ! apply boundary condition in all directions
        ! x -- periodic bc
        ! y -- free slip bc
        subroutine apply_parcel_bc(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            double precision                :: lower(2)
            double precision                :: upper(2)
            integer                         :: n

            lower = mesh%origin
            upper = mesh%origin + mesh%extent

            do n = 1, n_parcels
                ! horizontal bc
                call do_periodic(position, lower(1), upper(1), n)

                ! vertical bc
                call do_free_slip(position, velocity, lower(2), upper(2), n)
            enddo
        end subroutine apply_parcel_bc


        ! apply periodic boundary condition in x (horizontal) direction
        subroutine apply_periodic(position)
            double precision, intent(inout) :: position(:, :)
            double precision                :: lower, upper
            integer                         :: n

            lower = mesh%origin(1)
            upper = mesh%origin(1) + mesh%extent(1)

            do n = 1, n_parcels
                call do_periodic(position, lower, upper, n)
            enddo
        end subroutine apply_periodic


        ! apply periodic boundary condition in y (vertical) direction
        subroutine apply_free_slip(position, velocity)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            double precision                :: lower, upper
            integer                         :: n

            lower = mesh%origin(2)
            upper = mesh%origin(2) + mesh%extent(2)

            do n = 1, n_parcels
                call do_free_slip(position, velocity, lower, upper, n)
            enddo

        end subroutine apply_free_slip

        ! apply periodic bc on n-th parcel
        subroutine do_periodic(position, lower, upper, n)
            double precision, intent(inout) :: position(:, :)
            double precision, intent(in)    :: lower, upper
            integer,          intent(in)    :: n

            if (position(n, 1) > upper) then
                position(n, 1) = position(n, 1) - (upper - lower)
            endif

            if (position(n, 1) < lower) then
                position(n, 1) = position(n, 1) + (upper - lower)
            endif
        end subroutine do_periodic


        ! apply free slip bc on n-th parcel
        subroutine do_free_slip(position, velocity, lower, upper, n)
            double precision, intent(inout) :: position(:, :), velocity(:, :)
            double precision, intent(in)    :: lower, upper
            integer,          intent(in)    :: n

            if (position(n, 2) >= upper) then
                velocity(n, 2) = 0.0
                position(n, 2) = 2.0 * upper - position(n, 2)
            endif

            if (position(n, 2) <= lower) then
                velocity(n, 2) = 0.0
                position(n, 2) = 2.0 * lower - position(n, 2)
            endif
        end subroutine do_free_slip

end module parcel_bc
