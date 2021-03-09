module parcels
    implicit none

    double precision, allocatable, dimension(:) :: &
        x, y,        & ! positions
        dxdt, dydt,  & ! velocitues
        stretch

end module parcels
