! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : allow_larger_anisotropy, parcel, boundary
    use constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use mpi_communicator
    use mpi_layout, only : box, l_mpi_layout_initialised
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

    ! mesh spacing
    double precision, protected :: dx(3)

    ! inverse mesh spacing
    double precision, protected :: dxi(3)

    ! grid cell volume
    double precision, protected :: vcell

    ! inverse grid cell volume
    double precision, protected :: vcelli

    ! number of grid cells in each dimension
    integer :: nx, ny, nz

    ! total number of grid cells
    integer, protected :: ncell

    ! inverse of total number of grid cells
    double precision, protected :: ncelli

    ! total number of grid cells in horizontal plane (x, y)
    integer, protected :: nhcell

    ! inverse of total number of grid cells in horizontal plane (x, y)
    double precision, protected :: nhcelli

    ! total number of grid points
    integer, protected :: ngrid

    ! inverse of total number of grid points
    double precision, protected :: ngridi

    ! domain size
    double precision :: extent(3)

    ! inverse domain size
    double precision, protected :: extenti(3)

    ! domain volume
    double precision, protected :: vdomain

    ! inverse domain volume
    double precision, protected :: vdomaini

    ! domain centre
    double precision, protected :: center(3)

    ! domain half widths values
    double precision, protected :: hl(3)

    double precision, protected :: hli(3)

    ! domain origin
    double precision :: lower(3)

    ! domain upper boundary
    double precision, protected :: upper(3)

    ! minimum volume
    double precision, protected :: vmin

    ! upper bound for major semi-axis (used for splitting)
    double precision, protected :: amax

    ! maximum number of allowed parcels
    integer, protected :: max_num_parcels

    ! specifies if zeta is kept zero on a boundary;
    ! this also makes sure that dzeta/dt = 0 on a boundary
    logical, protected :: l_bndry_zeta_zero(2)

    contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        double precision :: msr

        if (.not. l_mpi_layout_initialised) then
            call mpi_print("MPI layout is not initialsed!")
            call mpi_stop
        endif

        upper = lower + extent

        extenti = one / extent

        dx = extent / dble((/nx, ny, nz/))
        dxi = one / dx;

        msr = maxval((/dxi(1) * dx(2), dxi(2) * dx(1),   &
                       dxi(1) * dx(3), dxi(3) * dx(1),   &
                       dxi(2) * dx(3), dxi(3) * dx(2)/))

        if (msr > two) then
            if (world%rank == world%root) then
                print *, '**********************************************************************'
                print *, '*                                                                    *'
                print *, '*   Warning: A mesh spacing ratio of more than 2 is not advisable!   *'
                print *, '*                                                                    *'
                print *, '**********************************************************************'
            endif

            if (.not. allow_larger_anisotropy) then
                call mpi_stop
            endif
        endif

        vdomain = product(extent)
        vdomaini = one / vdomain

        vcell = product(dx)
        vcelli = one / vcell

        nhcell = nx * ny
        nhcelli = one / dble(nhcell)

        ncell = nhcell * nz
        ncelli = one / dble(ncell)

        ! due to x periodicity it is only nx
        ngrid = nx * ny * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / parcel%min_vratio

        amax = minval(dx)

        max_num_parcels = int(box%halo_ncell * parcel%min_vratio * parcel%size_factor)

    end subroutine update_parameters


    subroutine set_zeta_boundary_flag(zeta)
        double precision, intent(in) :: zeta(-1:nz+1,                &
                                             box%hlo(2):box%hhi(2),  &
                                             box%hlo(1):box%hhi(1))
        double precision             :: rms_bndry(2), rms_interior, thres
        double precision             :: val(3)

        if (boundary%l_ignore_bndry_zeta_flag) then
            l_bndry_zeta_zero(:) = .false.
            if (world%rank == world%root) then
                print *, "WARNING: You allow the gridded vertical vorticity component"
                print *, "         at the boundaries to develop non-zero values."
                print *, "         Stop your simulation if this is not your intention."
            endif
        else
            ! rms interior
            val(1) = sum(zeta(1:nz-1, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)

            ! rms boundary
            val(2) = sum(zeta(0,      box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)
            val(3) = sum(zeta(nz,     box%lo(2):box%hi(2), box%lo(1):box%hi(1)) ** 2)

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               val(1:3),                &
                               3,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            rms_interior = dsqrt(val(1) * nhcelli / (nz-1))
            rms_bndry = dsqrt(val(2:3) * nhcelli)

            thres = boundary%zeta_tol * rms_interior + epsilon(rms_interior)

            l_bndry_zeta_zero(:) = (rms_bndry(:) < thres)
        endif

#if defined(ENABLE_VERBOSE) || !defined(NDEBUG)
        if ((.not. l_bndry_zeta_zero(1)) .and. (world%rank == world%root)) then
            print *, "WARNING: This simulation will not ensure that the gridded vertical"
            print *, "         vorticity component is zero at the lower boundary."
        endif

        if ((.not. l_bndry_zeta_zero(2)) .and. (world%rank == world%root)) then
            print *, "WARNING: This simulation will not ensure that the gridded vertical"
            print *, "         vorticity component is zero at the upper boundary."
        endif
#endif
    end subroutine set_zeta_boundary_flag


    subroutine read_zeta_boundary_flag(ncid)
        integer, intent(in)     :: ncid
        integer                 :: grp_ncid
        character(*), parameter :: name = 'parameters'

        if (boundary%l_ignore_bndry_zeta_flag) then
            l_bndry_zeta_zero(:) = .false.
            if (world%rank == world%root) then
                print *, "WARNING: You allow the gridded vertical vorticity component"
                print *, "         at the boundaries to develop non-zero values."
                print *, "         Stop your simulation if this is not your intention."
            endif
            return
        endif

        ncerr = nf90_inq_ncid(ncid, name, grp_ncid)

        if (ncerr == 0) then
            call read_netcdf_attribute(grp_ncid, 'l_lower_boundry_zeta_zero', l_bndry_zeta_zero(1))
            call read_netcdf_attribute(grp_ncid, 'l_upper_boundry_zeta_zero', l_bndry_zeta_zero(2))
        else
            if (world%rank == world%root) then
                print *, "WARNING: Could not find a '" // name // "' group in the provided"
                print *, "         NetCDF file."
                print *, "         Note this will result in the boundary zeta flags being"
                print *, "         set to false starting with parcels (which could well be"
                print *, "         undesirable)."
            endif
        endif

    end subroutine read_zeta_boundary_flag


    subroutine write_zeta_boundary_flag(ncid)
        integer, intent(in)     :: ncid
        integer                 :: grp_ncid
        character(*), parameter :: name = 'parameters'

        ncerr = nf90_def_grp(ncid, name, grp_ncid)

        call check_netcdf_error("Failed to create NetCDF group '" // name // "'.")

        call write_netcdf_attribute(grp_ncid, 'l_lower_boundry_zeta_zero', l_bndry_zeta_zero(1))
        call write_netcdf_attribute(grp_ncid, 'l_upper_boundry_zeta_zero', l_bndry_zeta_zero(2))

    end subroutine write_zeta_boundary_flag


    subroutine set_mesh_spacing(ext, nc)
        double precision, intent(in) :: ext(3)
        integer,          intent(in) :: nc(3)
        dx = ext / dble(nc)
    end subroutine set_mesh_spacing


    subroutine set_max_num_parcels(num)
        integer, intent(in) :: num
        max_num_parcels = num
    end subroutine set_max_num_parcels

#ifdef ENABLE_UNIT_TESTS
    subroutine set_vmin(val)
        double precision, intent(in) :: val
        vmin = val
    end subroutine set_vmin

    subroutine set_amax(val)
        double precision, intent(in) :: val
        amax = val
    end subroutine set_amax
#endif

end module parameters
