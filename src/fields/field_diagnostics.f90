! =============================================================================
!                   Field diagnostics that are written to HDF5.
! =============================================================================
module field_diagnostics
    use constants, only : zero
    use parameters, only : vcell, nx, nz, ngrid
    use fields
    use h5_utils
    use h5_writer
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


        subroutine write_h5_field_diagnostics(h5file_id)
            integer(hid_t), intent(in)    :: h5file_id
            double precision              :: rms_v, abserr_v
            integer                       :: max_npar, min_npar
#ifndef NDEBUG
            double precision              :: vol_sym_err
#endif

            !
            ! write diagnostics
            !
            rms_v = get_rms_volume_error()
            call write_h5_double_scalar_attrib(h5file_id, "rms volume error", rms_v)

            abserr_v = get_max_abs_normalised_volume_error()
            call write_h5_double_scalar_attrib(h5file_id, "max absolute normalised volume error", abserr_v)

            max_npar = maxval(nparg)
            call write_h5_int_scalar_attrib(h5file_id, "max num parcels per cell", max_npar)

            min_npar = minval(nparg)
            call write_h5_int_scalar_attrib(h5file_id, "min num parcels per cell", min_npar)

#ifndef NDEBUG
            vol_sym_err = maxval(dabs(sym_volg(0:nz, :)))
            call write_h5_double_scalar_attrib(h5file_id, "max symmetry volume error", &
                                               vol_sym_err)
#endif
        end subroutine write_h5_field_diagnostics

end module field_diagnostics
