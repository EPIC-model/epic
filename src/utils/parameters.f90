! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : box
    use constants
    implicit none

    ! mesh spacing
    double precision :: dx(2)

    ! inverse mesh spacing
    double precision :: dxi(2)

    ! grid cell volume, really area in 2D:
    double precision :: vcell

    ! number of grid cells in each dimension
    integer :: nx, nz

    ! total number of grid cells
    integer :: ncell

    ! total number of grid points
    integer :: ngrid

    ! domain size
    double precision :: extent(2)

    ! domain half widths values
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
        extent = dble(box%extent)
        lower  = dble(box%origin)

        dx = extent / dble(box%nc)
        dxi = one / dx;

        vcell = product(dx)

        nx = box%nc(1)
        nz = box%nc(2)

        ncell = nx * nz

        ! due to x periodicity it is only nx
        ngrid = nx * (nz + 1)

        ! domain
        hl = extent / two
        hli = one / hl

        ! we need to be symmetric around zero in x
        lower(1) = -hl(1)
        upper(1) = hl(1)

        lower(2) = dble(box%origin(2))
        upper(2) = dble(box%origin(2)) + extent(2)

    end subroutine update_parameters

end module parameters
