! =============================================================================
! This module contains modifiable physical parameters.
! =============================================================================
module phys_parameters
    use constants
    use phys_constants
    implicit none

    ![m] inverse condensation scale-height
    double precision, parameter :: lam_c = 0.01d0

    ![] see equation (5) of MPIC paper
    double precision, parameter :: glat = gravity * L_c / (c_p * theta_l0)

end module phys_parameters
