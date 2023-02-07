! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : allow_larger_anisotropy, parcel
    use constants
    implicit none

    ! mesh spacing
    double precision, protected :: dx(2)

    ! inverse mesh spacing
    double precision, protected :: dxi(2)

    ! grid cell volume, really area in 2D:
    double precision, protected :: vcell

    ! inverse grid cell volume
    double precision, protected :: vcelli

    ! number of grid cells in each dimension
    integer :: nx, nz

    ! total number of grid cells
    integer, protected :: ncell

    ! inverse of total number of grid cells
    double precision, protected :: ncelli

    ! total number of grid points
    integer, protected :: ngrid

    ! inverse of total number of grid points
    double precision, protected :: ngridi

    ! domain size
    double precision :: extent(2)

    ! inverse of domain size
    double precision, protected :: extenti(2)

    ! domain volume
    double precision, protected :: vdomain

    ! inverse domain volume
    double precision, protected :: vdomaini

    ! domain centre
    double precision, protected :: center(2)

    ! domain half widths values
    double precision, protected :: hl(2)

    double precision, protected :: hli(2)

    ! domain origin
    double precision :: lower(2)

    ! domain upper boundary
    double precision, protected :: upper(2)

    ! minimum volume
    double precision, protected :: vmin

    ! maximum volume
    double precision, protected :: vmax

    ! maximum number of allowed parcels
    integer, protected :: max_num_parcels

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters

        upper = lower + extent

        extenti = one / extent

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

        vdomain = product(extent)
        vdomaini = one / vdomain

        vcell = product(dx)
        vcelli = one / vcell

        ncell = nx * nz
        ncelli = one / dble(ncell)

        ! due to x periodicity it is only nx
        ngrid = nx * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / parcel%min_vratio
        vmax = vcell / parcel%max_vratio

        max_num_parcels = int(nx * nz * parcel%min_vratio * parcel%size_factor)

    end subroutine update_parameters

    subroutine set_mesh_spacing(ext, nc)
        double precision, intent(in) :: ext(2)
        integer,          intent(in) :: nc(2)
        dx = ext / dble(nc)
    end subroutine set_mesh_spacing

#ifdef ENABLE_UNIT_TESTS
    subroutine set_vmin(val)
        double precision, intent(in) :: val
        vmin = val
    end subroutine set_vmin


    subroutine set_vmax(val)
        double precision, intent(in) :: val
        vmax = val
    end subroutine set_vmax
#endif

end module parameters
