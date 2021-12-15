! =============================================================================
! This module contains non-modifiable physical constants.
! =============================================================================
module phys_constants
    use constants
    implicit none

    ![m/s**2]
    double precision, parameter :: gravity = 9.81d0

    ![J/kg] latent heat of condensation
    double precision, parameter :: L_v = 2.501e6

    ![J/(kg*K)] specific heat at constant pressure
    double precision, parameter :: c_p = 1005.0d0

    ![kg/m**3] saturation specific humidity at ground level
    double precision, parameter :: h_0 = 0.015d0

    ![K] mean liquid-water potential temperature
    double precision, parameter :: theta_l0 = 300.0d0

end module phys_constants
