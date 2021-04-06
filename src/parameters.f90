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

     ! grid cell volume, really area in 2D:
    double precision :: vcell

     ! number of grid cells
    integer :: ncell

    ! domain size
    double precision :: extent(2) = (/ pi, pi /)

    ! domain half widths and edge values:
    double precision :: hl(2)

    double precision :: hli(2)

    ! domain origin
    double precision :: lower(2)

    ! domain upper boundary
    double precision :: upper(2)

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        dx = extent / (grid - 1)
        dxi = one / dx;

        vcell = product(dx)

        ncell = product(grid - 1)

        ! domain
        hl = extent / two
        hli = one / hl
        lower = -hl
        upper = hl

    end subroutine update_parameters

end module parameters
