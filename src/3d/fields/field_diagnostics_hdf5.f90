! =============================================================================
!                   Write field diagnostics to HDF5.
! =============================================================================
module field_diagnostics_hdf5
    use field_diagnostics
    use h5_utils
    use h5_writer
    use h5_reader, only : get_num_steps
    use options, only : verbose, write_h5_options
    implicit none

    private

    integer :: n_writes

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname


    public :: create_h5_field_stats_file,   &
              write_h5_field_stats


    contains

        ! Create the field diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_h5_field_stats_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            h5fname =  basename // '_field_stats.hdf5'

            call exist_h5_file(h5fname, l_exist)

            if (l_restart .and. l_exist) then
                call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)
                call get_num_steps(h5file_id, n_writes)
                call close_h5_file(h5file_id)
                return
            endif

            n_writes = 0

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_scalar_attrib(h5file_id, 'output_type', 'field diagnostics')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id, lower, extent, (/nx, nz/))

            call close_h5_file(h5file_id)

        end subroutine create_h5_field_stats_file

        ! Write a step in the field diagnostic file.
        ! @param[in] t is the time
        subroutine write_h5_field_stats(t)
            double precision, intent(in)    :: t
            integer(hid_t)                  :: group
            character(:), allocatable       :: name
            logical                         :: created

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write field diagnostics to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            name = trim(get_step_group_name(n_writes))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif


            call write_h5_scalar_attrib(group, "t", t)

            !
            ! write diagnostics
            !
            call write_h5_scalar_attrib(group, "rms volume error", rms_v)

            call write_h5_scalar_attrib(group, "max absolute normalised volume error", abserr_v)

            call write_h5_scalar_attrib(group, "max num parcels per cell", max_npar)

            call write_h5_scalar_attrib(group, "min num parcels per cell", min_npar)

            call write_h5_scalar_attrib(group, "average num parcels per cell", avg_npar)

            call write_h5_scalar_attrib(group, "average num small parcels per cell", avg_nspar)

            ! increment counter
            n_writes = n_writes + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, n_writes)

            call close_h5_group(group)

            call close_h5_file(h5file_id)

        end subroutine write_h5_field_stats

end module field_diagnostics_hdf5
