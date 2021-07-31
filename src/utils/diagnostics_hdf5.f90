module diagnostics_hdf5
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics, only : write_h5_parcel_diagnostics
    use field_diagnostics, only : write_h5_field_diagnostics
    use options, only : verbose
    use hdf5
    use h5_utils
    use h5_writer
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: hdf5_diagnostics_timer

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    private :: h5file_id, h5fname

    contains

        subroutine create_h5_diagnostics_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical                  :: overwrite

            h5fname =  basename // '_diagnostics.hdf5'

            call create_h5_file(h5fname, h5file_id)

            call write_h5_char_scalar_attrib(h5file_id, 'output_type', 'diagnostics')

            call close_h5_file(h5file_id)

        end subroutine create_h5_diagnostics_file

        subroutine write_h5_diagnostics_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            integer(hid_t)                  :: group, pgroup, fgroup
            character(:), allocatable       :: name
            logical                         :: created
            character(len=6)                :: pname = 'parcel'
            character(len=5)                :: fname = 'field'

            call start_timer(hdf5_diagnostics_timer)


#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write diagnostics to h5"
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
            ! write parcel diagnostics
            !
            call create_h5_group(group, pname, pgroup, created)

            if (.not. created) then
                call open_h5_group(group, pname, pgroup)
            endif

            call write_h5_parcel_diagnostics(pgroup)

            call close_h5_group(pgroup)

            !
            ! write field diagnostics
            !
            call create_h5_group(group, fname, fgroup, created)

            if (.not. created) then
                call open_h5_group(group, fname, fgroup)
            endif

            call write_h5_field_diagnostics(fgroup)

            call close_h5_group(fgroup)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_group(group)

            call close_h5_file(h5file_id)

            call stop_timer(hdf5_diagnostics_timer)

        end subroutine write_h5_diagnostics_step

end module diagnostics_hdf5
