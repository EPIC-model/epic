! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : allow_larger_anisotropy
    use constants
    use h5_writer
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
        hl = extent / two
        hli = one / hl

    end subroutine update_parameters


    subroutine write_h5_parameters(fname)
        character(*), intent(in) :: fname
!         integer(hid_t)           :: group
!
!         call open_h5_file(fname)
!         group = open_h5_group("box")
!         call write_h5_double_vector_attrib(group, "extent", extent)
!         call write_h5_double_vector_attrib(group, "origin", lower)
!         call write_h5_integer_vector_attrib(group, "ncells", (/nx, nz/))
!         call h5gclose_f(group, h5err)
!         call close_h5_file
    end subroutine write_h5_parameters

end module parameters
