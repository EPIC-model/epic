! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use macros, only : EPIC_VECTOR
    use options, only : allow_larger_anisotropy, parcel
    use constants
    implicit none

    ! mesh spacing
    double precision :: dx(ndim)

    ! inverse mesh spacing
    double precision :: dxi(ndim)

    ! grid cell volume, really area in 2D:
    double precision :: vcell

    ! inverse grid cell volume
    double precision :: vcelli

    ! number of grid cells in each dimension
    integer :: nx, ny, nz

    ! total number of grid cells
    integer :: ncell

    ! inverse of total number of grid cells
    double precision :: ncelli

    ! total number of grid points
    integer :: ngrid

    ! inverse of total number of grid points
    double precision :: ngridi

    ! domain size
    double precision :: extent(ndim)

    ! domain centre
    double precision :: center(ndim)

    ! domain half widths values
    double precision :: hl(ndim)

    double precision :: hli(ndim)

    ! domain origin
    double precision :: lower(ndim)

    ! domain upper boundary
    double precision :: upper(ndim)

    ! minimum volume
    double precision :: vmin

    ! maximum volume
    double precision :: vmax

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        double precision :: msr ! mesh spacing ratio

        upper = lower + extent

        dx = extent / dble(EPIC_VECTOR(nx, ny, nz))
        dxi = one / dx;

#ifdef ENABLE_3D
        msr = maxval((/dxi(1) * dx(2), dxi(2) * dx(1),   &
                       dxi(1) * dx(3), dxi(3) * dx(1),   &
                       dxi(2) * dx(3), dxi(3) * dx(2)/))
#else
        msr = max(dxi(1) * dx(2), dxi(2) * dx(1))
#endif

        if (msr > two) then
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
        vcelli = one / vcell

        ncell = nx * ny * nz
        ncelli = one / dble(ncell)

        ! due to x periodicity it is only nx (and also ny in 3D)
        ngrid = nx * ny * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / parcel%min_vratio
        vmax = vcell / parcel%max_vratio

    end subroutine update_parameters
end module parameters
