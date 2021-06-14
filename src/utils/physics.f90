! =============================================================================
! This module contains non-modifiable physical constants.
! =============================================================================
module physics
    implicit none

    double precision, parameter :: gravity = 9.81d0   ! m / s**2

    ! surface saturation humidity
    double precision, parameter :: q0 = 0.0d0 ! FIXME

    ! see equation (5) of MPIC paper
    double precision, parameter :: glat = 0.0d0 ! FIXME

end module physics
