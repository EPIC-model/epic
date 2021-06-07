! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero
    use parcel_container, only : parcels, n_parcels
    use hdf5
    use writer, only : h5file,                         &
                       h5err,                          &
                       write_h5_double_scalar_attrib,  &
                       open_h5_group,                  &
                       get_step_group_name
    implicit none

    private

    ! peref : potential energy reference
    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: peref, pe, ke

    public :: init_parcel_diagnostics, &
              write_h5_parcel_diagnostics

    contains

        ! compute the reference potential energy
        subroutine init_parcel_diagnostics

            !TODO
            peref = zero

        end subroutine init_parcel_diagnostics


        subroutine calculate_diagnostics
            integer          :: n
            double precision :: b, z, v, vel(2), vol

            ! reset
            ke = zero
            pe = zero

            do n = 1, n_parcels

                vel = parcels%velocity(n, :)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(n, 2)

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                ! potential energy
                pe = pe - b * z * vol
            enddo

            ke = 0.5d0 * ke
            pe = pe - peref
        end subroutine calculate_diagnostics


        subroutine write_h5_parcel_diagnostics(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            integer(hid_t)                :: step_group
            character(:), allocatable     :: step
            character(:), allocatable     :: name

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/diagnostics"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            call calculate_diagnostics

            !
            ! write diagnostics
            !
            call write_h5_double_scalar_attrib(group, "potential energy", pe)
            call write_h5_double_scalar_attrib(group, "kinetic energy", ke)
            call write_h5_double_scalar_attrib(group, "total energy", ke + pe)


            ! close all
            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)
        end subroutine write_h5_parcel_diagnostics



end module parcel_diagnostics
