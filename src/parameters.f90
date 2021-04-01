! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : grid
    use constants
    implicit none

    ! mesh spacing
    double precision :: dx(2)

    ! inverse mesh spacing
    double precision :: dxi(2)

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        dx = extent / (grid - 1)
        dxi = one / dx;
    end subroutine update_parameters

end module parameters
