! =============================================================================
!                   Field diagnostics that are written to HDF5.
! =============================================================================
module field_diagnostics
    use constants, only : zero
    use parameters, only : vcell, nx, nz, ngrid, ncell
    use fields
    use h5_utils
    use h5_writer
    use timer, only : start_timer, stop_timer
    implicit none

    private

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    integer :: hdf5_field_stat_timer

    public :: create_h5_field_stats_file, &
              write_h5_field_stats_step,  &
              hdf5_field_stat_timer


    contains

        subroutine create_h5_field_stats_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical,      intent(in) :: overwrite

            h5fname =  basename // '_field_stats.hdf5'

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_char_scalar_attrib(h5file_id, 'output_type', 'field diagnostics')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id)

            call close_h5_file(h5file_id)

        end subroutine create_h5_field_stats_file

        function get_max_abs_normalised_volume_error() result(err)
            double precision :: err
            err = maxval(abs(volg(0:nz, :)  - vcell)) / vcell
        end function get_max_abs_normalised_volume_error

        function get_rms_volume_error() result(rms)
            double precision :: rms ! rms volume error
            double precision :: sqerrsum

            ! do not take halo cells into account
            sqerrsum = sum((volg(0:nz, :) - vcell) ** 2)

            rms = dsqrt(sqerrsum / dble(ngrid)) / vcell
        end function get_rms_volume_error


        subroutine write_h5_field_stats_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            integer(hid_t)                  :: group
            character(:), allocatable       :: name
            logical                         :: created
            double precision                :: rms_v, abserr_v, res
            integer                         :: max_npar, min_npar
#ifndef NDEBUG
            double precision                :: vol_sym_err
#endif

            call start_timer(hdf5_field_stat_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write field diagnostics to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif


            call write_h5_double_scalar_attrib(group, "t", t)

            call write_h5_double_scalar_attrib(group, "dt", dt)

            !
            ! write diagnostics
            !
            rms_v = get_rms_volume_error()
            call write_h5_double_scalar_attrib(group, "rms volume error", rms_v)

            abserr_v = get_max_abs_normalised_volume_error()
            call write_h5_double_scalar_attrib(group, "max absolute normalised volume error", abserr_v)

            max_npar = maxval(nparg(0:nz-1, :))
            call write_h5_int_scalar_attrib(group, "max num parcels per cell", max_npar)

            min_npar = minval(nparg(0:nz-1, :))
            call write_h5_int_scalar_attrib(group, "min num parcels per cell", min_npar)

            res = sum(nparg(0:nz-1, :)) / dble(ncell)
            call write_h5_double_scalar_attrib(group, "average num parcels per cell", res)

            res = sum(nsparg(0:nz-1, :)) / dble(ncell)
            call write_h5_double_scalar_attrib(group, "average num small parcels per cell", res)

#ifndef NDEBUG
            vol_sym_err = maxval(dabs(sym_volg(0:nz, :)))
            call write_h5_double_scalar_attrib(group, "max symmetry volume error", &
                                               vol_sym_err)
#endif

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_group(group)

            call close_h5_file(h5file_id)

            call stop_timer(hdf5_field_stat_timer)

        end subroutine write_h5_field_stats_step

end module field_diagnostics
