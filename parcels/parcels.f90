module parcels
    implicit none

    double precision, allocatable, dimension(:) :: &
        x, y,        & ! positions
        dxdt, dydt,  & ! velocitues
        stretch

    contains

        subroutine split(threshold)
            double precision, intent(in) :: threshold


        end subroutine split

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcels
