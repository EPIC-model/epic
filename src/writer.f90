!
!
! Fortran HDF5 examples at
! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
module writer
    use hdf5
    implicit none

    ! h5 file handle
    integer(HID_T) :: h5file = 0

    ! if non-zero an error occurred
    integer :: h5err = 0

    contains
        subroutine open_h5_file(filename)
            character(*), intent(in) :: filename

            if (h5file .eq. 0) then
                call h5fcreate_f(filename, H5F_ACC_TRUNC_F, h5file, h5err)
            endif

            call h5open_f(h5err)

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

        subroutine write_h5_dataset_1d(group, name, data)
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:1) :: dims

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(1, dims, dataspace, h5err)

            ! create the dataset
            call h5dcreate_f(h5file, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)

            ! close all
            call h5dclose_f(dset , h5err)
            call h5sclose_f(dataspace, h5err)
        end subroutine write_h5_dataset_1d

        subroutine write_h5_dataset_2d(group, name, data)
            ! 12 March 2021
            ! https://stackoverflow.com/questions/48816383/passing-character-strings-of-different-lengths-to-functions-in-fortran
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:, :)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:2) :: dims

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(2, dims, dataspace, h5err)

            ! create the dataset
            call h5dcreate_f(h5file, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)

            ! close all
            call h5dclose_f(dset , h5err)
            call h5sclose_f(dataspace, h5err)
        end subroutine write_h5_dataset_2d


        subroutine write_h5_scalar_attrib(iter, name, val)
            integer,          intent(in)     :: iter ! iteration
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: val
            integer(hid_t)                   :: group
            character(:), allocatable        :: grn
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1

            grn = trim(get_step_group_name(iter))

            ! create or open group
            group = open_h5_group(grn)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_DOUBLE, space, attr, h5err)

            call h5awrite_f(attr, H5T_NATIVE_DOUBLE, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
            call h5gclose_f(group, h5err)
        end subroutine write_h5_scalar_attrib


        ! open existing group or create one
        function open_h5_group(name) result(group)
            character(*),   intent(in) :: name
            integer(hid_t)             :: group
            logical                    :: link_exists = .false.

            call h5lexists_f(h5file, name, link_exists, h5err)

            if (link_exists) then
                call h5gopen_f(h5file, name, group, h5err)
            else
                call h5gcreate_f(h5file, name, group, h5err)
            endif
        end function open_h5_group

        subroutine close_h5_group(group)
            integer(hid_t), intent(in) :: group
            call h5gclose_f(group, h5err)
        end subroutine close_h5_group


        ! convert iteration number to string
        function get_step_group_name(iter) result(name)
            integer, intent(in) :: iter
            ! 12 March 2021
            ! https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
            character(len=32) :: name

            write(name, fmt='(I10.10)') iter
            name = 'step#' // trim(name)
        end function get_step_group_name

end module writer
