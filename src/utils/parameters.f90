! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : allow_larger_anisotropy, parcel
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

    ! domain centre
    double precision :: center(2)

    ! domain half widths values
    double precision :: hl(2)

    double precision :: hli(2)

    ! domain origin
    double precision :: lower(2)

    ! domain upper boundary
    double precision :: upper(2)

    ! minimum volume
    double precision :: vmin

    ! maximum volume
    double precision :: vmax

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters

        upper = lower + extent

        dx = extent / dble((/nx, nz/))
        dxi = one / dx;

        if (max(dxi(1) * dx(2), dxi(2) * dx(1)) > two) then
            print *, '**********************************************************************'
            print *, '*                                                                    *'
            print *, '*   Warning: A mesh spacing ratio of more than 2 is not advisable!   *'
            print *, '*                                                                    *'
            print *, '**********************************************************************'

            if (.not. allow_larger_anisotropy) then
                stop
            endif
        endif

        vcell = product(dx)

        ncell = nx * nz

        ! due to x periodicity it is only nx
        ngrid = nx * (nz + 1)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / dble(parcel%vmin_fraction)
        vmax = vcell / dble(parcel%vmax_fraction)

    end subroutine update_parameters
end module parameters
