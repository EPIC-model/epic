! =============================================================================
! This module contains modifiable physical parameters.
! =============================================================================
module phys_parameters
    use constants
    use options, only : l_coriolis, lat_degrees, ang_vel
    use phys_constants
    implicit none

    ![m] inverse condensation scale-height
    double precision, parameter :: lam_c = 0.001d0

    ![] see equation (5) of MPIC paper
    double precision, parameter :: glat = gravity * L_v / (c_p * theta_l0)

    double precision, parameter :: glati = one / glat

    !FIXME comment
    double precision :: lat_ref

    ! Coriolis frequency
    double precision :: f_cor

    ! component of the planetary vorticity in the y direction
    double precision :: ft_cor

    contains

        subroutine update_phys_parameters
            if (l_coriolis) then
                lat_ref = lat_degrees * deg2rad
                f_cor  = two * ang_vel * dsin(lat_ref)
                ft_cor = two * ang_vel * dcos(lat_ref)
            else
                f_cor  = zero
                ft_cor = zero
            endif
        end subroutine update_phys_parameters

end module phys_parameters
