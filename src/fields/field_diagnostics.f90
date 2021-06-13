! =============================================================================
!                   Field diagnostics that are written to HDF5.
! =============================================================================
module field_diagnostics
    use constants, only : zero
    use parameters, only : vcell, nx, nz, ngrid
    use fields
!     use taylorgreen, only : get_flow_vorticity
    use hdf5
    use writer, only : h5file,                         &
                       h5err,                          &
                       write_h5_double_scalar_attrib,  &
                       write_h5_integer_scalar_attrib, &
                       open_h5_group,                  &
                       get_step_group_name
    implicit none


    contains

        function get_max_abs_normalised_volume_error() result(err)
            double precision :: err
            err = maxval(abs(volg(0:nz, 0:nx-1)  - vcell)) / vcell
        end function get_max_abs_normalised_volume_error

        function get_rms_volume_error() result(rms)
            double precision :: rms ! rms volume error
            double precision :: sqerrsum

            ! do not take halo cells into account
            ! x indices run from 0 to nx-1
            sqerrsum = sum((volg(0:nz, 0:nx-1) - vcell) ** 2)

            rms = dsqrt(sqerrsum / dble(ngrid)) / vcell
        end function get_rms_volume_error

!         function get_rms_vorticity_error() result(rms)
!             double precision :: rms
!             double precision :: sqerrsum, vexact, pos(2)
!             integer          :: i, j
!
!             sqerrsum = zero
!             do j = 0, nz
!                 do i = 0, nx-1
!                     call get_position(i, j, pos)
!                     vexact = get_flow_vorticity(pos)
!                     sqerrsum = sqerrsum + (vortg(j, i) - vexact) ** 2
!                 enddo
!             enddo
!
!             rms = dsqrt(sqerrsum / dble(ngrid))
!
!         end function get_rms_vorticity_error


        subroutine write_h5_field_diagnostics(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            integer(hid_t)                :: step_group
            character(:), allocatable     :: step
            character(:), allocatable     :: name
            double precision              :: rms_v, abserr_v
            integer                       :: max_npar, min_npar

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

            abserr_v = get_max_abs_normalised_volume_error()
            call write_h5_double_scalar_attrib(group, "max absolute normalised volume error", abserr_v)

            max_npar = maxval(nparg)
            call write_h5_integer_scalar_attrib(group, "max num parcels per cell", max_npar)

            min_npar = minval(nparg)
            call write_h5_integer_scalar_attrib(group, "min num parcels per cell", min_npar)

!             rms_v = get_rms_vorticity_error()
!             call write_h5_double_scalar_attrib(group, "rms vorticity error", rms_v)


            ! close all
            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)
        end subroutine write_h5_field_diagnostics

end module field_diagnostics
