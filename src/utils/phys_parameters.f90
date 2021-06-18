! =============================================================================
! This module contains modifiable physical parameters.
! =============================================================================
module phys_parameters
    use phys_constants
    implicit none

    ![g/m**3] surface saturation humidity
    double precision, parameter :: q0_c = 0.0d0 ! FIXME

    ![J/g] latent heat of condensation
    double precision, parameter :: L_c = 0.0d0 ! FIXME

    ![J/(g*K)] specific heat at constant pressure
    double precision, parameter :: cp_c = 1.0d0 ! FIXME

    ![m] inverse condensation scale-height
    double precision, parameter :: lam_c = 0.0d0 ! FIXME

    ![K] liquid-water potential temperature reference
    double precision, parameter :: theta_c = 1.0d0 ! FIXME

    ![] see equation (5) of MPIC paper
    double precision, parameter :: glat_c = g_c * L_c / (cp_c * theta_c)

end module phys_parameters
