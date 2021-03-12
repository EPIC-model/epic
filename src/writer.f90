!
!
! Fortran HDF5 examples at
! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
module writer
    use hdf5
    implicit none

    !character(:), allocatable :: filename

    ! h5 file handle
    integer(HID_T) :: h5file

    ! if non-zero an error occurred
    integer :: h5err = 0

    contains
        subroutine h5_create_file
            call h5fcreate_f('test.h5', H5F_ACC_TRUNC_F, h5file, h5err)

            if (h5err .ne. 0) then
                print *, "Could not create the H5 file."
                stop
            endif
        end subroutine h5_create_file


        subroutine h5_open_file
            call h5open_f(h5err)

            if (h5err .ne. 0) then
                print *, "Opening the H5 file failed."
                stop
            endif
        end subroutine h5_open_file

        subroutine h5_close_file
            call h5close_f(h5err)

            if (h5err .ne. 0) then
                print *, "Closing the H5 file failed."
                stop
            endif
        end subroutine h5_close_file

        subroutine write_h5_dataset_1d(name, data)
            character(*), intent(in) :: name
            double precision, intent(in) :: data(:)
            integer(hid_t)  :: dset, dataspace
            integer(hsize_t), dimension(1:1) :: dims

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(1, dims, dataspace, h5err)

            ! create the dataset
            call h5dcreate_f(h5file, name, H5T_NATIVE_DOUBLE, dataspace, dset, h5err)

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)

            ! close all
            call h5dclose_f(dset , h5err)
            call h5sclose_f(dataspace, h5err)
        end subroutine

        subroutine write_h5_dataset_2d(name, data)
            ! 12 March 2021
            ! https://stackoverflow.com/questions/48816383/passing-character-strings-of-different-lengths-to-functions-in-fortran
            character(*), intent(in) :: name
            double precision, intent(in) :: data(:, :)
            integer(hid_t)  :: dset, dataspace
            integer(hsize_t), dimension(1:2) :: dims

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(2, dims, dataspace, h5err)

            ! create the dataset
            call h5dcreate_f(h5file, name, H5T_NATIVE_DOUBLE, dataspace, dset, h5err)

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)

            ! close all
            call h5dclose_f(dset , h5err)
            call h5sclose_f(dataspace, h5err)
        end subroutine


end module writer
