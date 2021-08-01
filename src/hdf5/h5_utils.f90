! =============================================================================
!               HDF5 utilities common to h5_reader and h5_writer
! =============================================================================
module h5_utils
    use hdf5
    implicit none

    ! if non-zero an error occurred
    integer :: h5err = 0

    logical :: is_initialised = .false.

    private :: is_initialised

    contains

        subroutine initialise_hdf5
            call h5open_f(h5err)

            call check_h5_error("Failed to start hdf5 properly.")

            is_initialised = .true.

        end subroutine initialise_hdf5


        subroutine finalise_hdf5
            if (.not. is_initialised) then
                return
            endif

            call h5close_f(h5err)

            call check_h5_error("Failed to terminated hdf5 properly.")
        end subroutine finalise_hdf5


        subroutine create_h5_file(h5fname, overwrite, h5file_id)
            character(*),   intent(in)  :: h5fname
            logical,        intent(in)  :: overwrite
            integer(hid_t), intent(out) :: h5file_id
            logical                     :: exists = .true.

            ! check whether file exists
            inquire(file=h5fname, exist=exists)

            if (exists .and. overwrite) then
                call delete_h5_file(trim(h5fname))
            else if (exists) then
                print *, "File '" // trim(h5fname) // "' already exists. Exiting."
                stop
            endif

            call h5fcreate_f(h5fname, H5F_ACC_TRUNC_F, h5file_id, h5err)
            call check_h5_error("Failed to create hdf5 file'" // trim(h5fname) // "'.")
        end subroutine create_h5_file

        subroutine delete_h5_file(h5fname)
            character(*), intent(in) :: h5fname
            integer                  :: stat

            ! 15 June 2021
            ! https://stackoverflow.com/questions/18668832/how-delete-file-from-fortran-code
            open(unit=1234, iostat=stat, file=h5fname, status='old')
            if (stat == 0) then
                close(1234, status='delete')
            endif
        end subroutine delete_h5_file


        subroutine open_h5_file(h5fname, access_flag, h5file_id)
            character(*),   intent(in)  :: h5fname
            integer,        intent(in)  :: access_flag ! H5F_ACC_RDWR_F or ! H5F_ACC_RDONLY_F
            integer(hid_t), intent(out) :: h5file_id

            call h5fopen_f(h5fname, access_flag, h5file_id, h5err)

            call check_h5_error("Opening the hdf5 file failed.")
        end subroutine open_h5_file


        subroutine close_h5_file(h5file_id)
            integer(hid_t), intent(in) :: h5file_id

            call h5fclose_f(h5file_id, h5err)

            call check_h5_error("Closing the hdf5 file failed.")
        end subroutine close_h5_file


        subroutine create_h5_group(h5file_id, name, group, created)
            integer(hid_t),            intent(in)  :: h5file_id
            character(*),              intent(in)  :: name
            integer(hid_t),            intent(out) :: group
            logical,        optional,  intent(out) :: created
            logical                                :: link_exists = .false.

            call h5lexists_f(h5file_id, name, link_exists, h5err)

            call check_h5_error("Failed to check if hdf5 link exists.")

            if (link_exists) then
                if (present(created)) then
                    created = .false.
                endif
                return
            endif

            if (present(created)) then
                created = .true.
            endif

            call h5gcreate_f(h5file_id, name, group, h5err)
            call check_h5_error("Failed to create hdf5 group.")
        end subroutine create_h5_group


        ! open existing group
        subroutine open_h5_group(h5file_id, name, group)
            integer(hid_t), intent(in)  :: h5file_id
            character(*),   intent(in)  :: name
            integer(hid_t), intent(out) :: group
            logical                     :: link_exists = .false.

            call h5lexists_f(h5file_id, name, link_exists, h5err)

            call check_h5_error("Failed to check if hdf5 link exists.")

            if (link_exists) then
                call h5gopen_f(h5file_id, name, group, h5err)
                call check_h5_error("Failed to open hdf5 group.")
            else
                print *, "Group '", name, "' does not exist!"
                stop
            endif
        end subroutine open_h5_group

        subroutine close_h5_group(group)
            integer(hid_t) :: group
            call h5gclose_f(group, h5err)
            call check_h5_error("Failed to close hdf5 group.")
        end subroutine close_h5_group


        subroutine check_h5_error(msg)
            character(*), intent(in) :: msg
#ifndef NDEBUG
            if (h5err .ne. 0) then
                print *, msg
                stop
            endif
#endif
        end subroutine check_h5_error
end module h5_utils
