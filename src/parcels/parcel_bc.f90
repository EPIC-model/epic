module parcel_bc
    use parameters, only : mesh
    use parcel_container, only : n_parcels
    implicit none

    private :: apply_free_slip,    &
               apply_periodic

    contains

        subroutine apply_parcel_bc(position, velocity)
            double precision, intent(inout) :: position(:, :)
            double precision, intent(inout) :: velocity(:, :)
            double precision                :: lower(2)
            double precision                :: upper(2)
            integer                         :: i

            lower = mesh%origin
            upper = mesh%origin + mesh%extent

            do i = 1, 2
                if (mesh%bc(i) == "free slip") then
                    call apply_free_slip(position(:, i), &
                                         velocity(:, i), &
                                         lower(i), upper(i))
                else if (mesh%bc(i) == "periodic") then
                    call apply_periodic(position(:, i), lower(i), upper(i))
                else
                    print *, "No boundary condition named '", mesh%bc(i), "'."
                    stop
                endif
            enddo
        end subroutine apply_parcel_bc


        subroutine apply_periodic(position, lower, upper)
            double precision, intent(inout) :: position(:)
            double precision, intent(in)    :: lower, upper
            integer                         :: n

            do n = 1, n_parcels
                if (position(n) > upper) then
                    position(n) = position(n) - (upper - lower)
                endif

                if (position(n) < lower) then
                    position(n) = position(n) + (upper - lower)
                endif
            enddo

        end subroutine apply_periodic


        subroutine apply_free_slip(position, velocity, lower, upper)
            double precision, intent(inout) :: position(:), velocity(:)
            double precision, intent(in)    :: lower, upper
            integer                         :: n

            do n = 1, n_parcels
                if (position(n) >= upper) then
                    velocity(n) = 0.0
                    position(n) = 2.0 * upper - position(n)
                endif

                if (position(n) <= lower) then
                    velocity(n) = 0.0
                    position(n) = 2.0 * lower - position(n)
                endif
            enddo

        end subroutine apply_free_slip

end module parcel_bc
