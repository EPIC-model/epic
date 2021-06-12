! =============================================================================
!                               HDF5 reader
!   This module reads in datatsets. For further information about the
!   HDF5 API please visit
!   https://support.hdfgroup.org/HDF5/doc1.6/RM_H5Front.html
! =============================================================================
module reader
    use hdf5
    implicit none

    ! h5 file handle
    integer(hid_t) :: h5file = 0

    ! if non-zero an error occurred
    integer :: h5err = 0

    private :: h5file, h5err

    contains

        subroutine open_h5_file(filename)
            character(*), intent(in) :: filename

            call h5open_f(h5err)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, h5file, h5err)

            if (h5err .ne. 0) then
                print *, "Opening the H5 file failed."
                stop
            endif
        end subroutine open_h5_file

        subroutine close_h5_file
            call h5close_f(h5err)

            if (h5err .ne. 0) then
                print *, "Closing the H5 file failed."
                stop
            endif
        end subroutine close_h5_file


        function has_dataset(name) result(link_exists)
            character(*), intent(in) :: name
            logical                  :: link_exists
            link_exists = .false.
            call h5lexists_f(h5file, name, link_exists, h5err)
        end function has_dataset

        subroutine read_h5_dataset_2d(name, buffer)
            character(*),                  intent(in)  :: name
            double precision, allocatable, intent(out) :: buffer(:, :)
            integer(hid_t)                             :: dataspace_id   ! Dataspace identifier
            integer(hid_t)                             :: dset_id
            integer                                    :: rank
            integer(hsize_t)                           :: dims(2), maxdims(2)

            call h5dopen_f(h5file, name, dset_id, h5err)
            call h5dget_space_f(dset_id, dataspace_id, h5err)

            call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)

            call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, h5err)

            if (.not. sum(dims -  maxdims) == 0) then
                print *, "Dimensions do not agree."
                stop
            endif

            allocate(buffer(dims(1), dims(2)))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)

            call h5dclose_f(dset_id, h5err)
        end subroutine read_h5_dataset_2d

end module reader
