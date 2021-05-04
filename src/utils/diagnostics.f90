! =============================================================================
!           This module contains diagnostics that are written to HDF5.
! =============================================================================
module diagnostics
    use parcel_container
    use parameters, only : vcell, ncell, nx, nz
    use fields
    use hdf5
    use writer, only : h5file,                        &
                       h5err,                         &
                       write_h5_double_scalar_attrib, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none


    contains

        function get_max_abs_volume_error() result(err)
            double precision :: err
            err = maxval(abs(volg(0:nz, 0:nx-1, :)  - vcell))
        end function get_max_abs_volume_error

        function get_rms_volume_error() result(rms)
            double precision :: rms ! rms volume error
            double precision :: sqerrsum

            ! do not take halo cells into account
            ! x indices run from 0 to nx-1
            sqerrsum = sum((volg(0:nz, 0:nx-1, :) - vcell) ** 2)

            rms = sqrt(sqerrsum / dble(ncell)) / vcell
        end function get_rms_volume_error

        subroutine write_h5_diagnostics(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            integer(hid_t)                :: step_group
            character(:), allocatable     :: step
            character(:), allocatable     :: name
            double precision              :: rms_v, abserr_v

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/diagnostics"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            !
            ! write diagnostics
            !

            rms_v = get_rms_volume_error()
            call write_h5_double_scalar_attrib(group, "rms volume error", rms_v)

            abserr_v = get_max_abs_volume_error()
            call write_h5_double_scalar_attrib(group, "max absolute volume error", abserr_v)


            ! close all
            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)
        end subroutine write_h5_diagnostics

end module diagnostics
