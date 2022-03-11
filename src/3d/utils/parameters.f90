! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : allow_larger_anisotropy, parcel
    use constants
    implicit none

    ! mesh spacing
    double precision :: dx(3)

    ! inverse mesh spacing
    double precision :: dxi(3)

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
    double precision :: extent(3)

    ! domain centre
    double precision :: center(3)

    ! domain half widths values
    double precision :: hl(3)

    double precision :: hli(3)

    ! domain origin
    double precision :: lower(3)

    ! domain upper boundary
    double precision :: upper(3)

    ! minimum volume
    double precision :: vmin

    ! maximum volume
    double precision :: vmax

    ! maximum number of allowed parcels
    integer :: max_num_parcels

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        double precision :: msr

        upper = lower + extent

        dx = extent / dble((/nx, ny, nz/))
        dxi = one / dx;

        msr = maxval((/dxi(1) * dx(2), dxi(2) * dx(1),   &
                       dxi(1) * dx(3), dxi(3) * dx(1),   &
                       dxi(2) * dx(3), dxi(3) * dx(2)/))

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

        ! due to x periodicity it is only nx
        ngrid = nx * ny * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / parcel%min_vratio
        vmax = vcell / parcel%max_vratio

        max_num_parcels = nx * ny * nz * parcel%n_per_cell * parcel%size_factor

    end subroutine update_parameters
end module parameters
