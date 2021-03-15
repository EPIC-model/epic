!
!
! Fortran HDF5 examples at
! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
module writer
    use hdf5
    implicit none

    ! h5 file handle
    integer(hid_t) :: h5file = 0

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

            if (size(data) == 1) then
                print *, "Error in 'write_h5_dataset_1d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

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

            if (size(data) == 1) then
                print *, "Error in 'write_h5_dataset_2d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

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

        subroutine write_h5_dataset_3d(group, name, data)
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:, :, :)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:3) :: dims

            if (size(data) == 1) then
                print *, "Error in 'write_h5_dataset_3d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(3, dims, dataspace, h5err)

            ! create the dataset
            call h5dcreate_f(h5file, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)

            ! close all
            call h5dclose_f(dset , h5err)
            call h5sclose_f(dataspace, h5err)
        end subroutine write_h5_dataset_3d


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

        ! convert iteration number to string
        function get_step_group_name(iter) result(name)
            integer, intent(in) :: iter
            ! 12 March 2021
            ! https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
            character(len=32) :: name

            write(name, fmt='(I10.10)') iter
            name = 'step#' // trim(name)
        end function get_step_group_name


        subroutine write_h5_double_vector_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: val(:)
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims

            dims = size(val)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_DOUBLE, space, attr, h5err)

            call h5awrite_f(attr, H5T_NATIVE_DOUBLE, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
        end subroutine write_h5_double_vector_attrib

        subroutine write_h5_double_scalar_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_DOUBLE, space, attr, h5err)

            call h5awrite_f(attr, H5T_NATIVE_DOUBLE, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
        end subroutine write_h5_double_scalar_attrib

        subroutine write_h5_integer_scalar_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
        end subroutine write_h5_integer_scalar_attrib


        subroutine write_h5_integer_vector_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: val(:)
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims

            dims = size(val)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
        end subroutine write_h5_integer_vector_attrib


        subroutine write_h5_character_scalar_attrib(group, name, val)
            integer(hid_t),    intent(in)     :: group
            character(*),      intent(in)     :: name
            character(len=16), intent(in)     :: val
            integer(hid_t)                    :: attr, space, filetype, memtype
            integer(hsize_t),  dimension(1:1) :: dims = 1
            integer(size_t),   parameter      :: sdim = 16

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! type in file
            call h5tcopy_f(H5T_C_S1, filetype, h5err)
            call h5tset_size_f(filetype, sdim+1, h5err)

            ! type in run
            call h5tcopy_f(H5T_FORTRAN_S1, memtype, h5err)
            call h5tset_size_f(memtype, sdim, h5err)

            ! create the dataset
            call h5acreate_f(group, name, filetype, space, attr, h5err)

            ! write attribte
            call h5awrite_f(attr, memtype, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
            call h5tclose_f(filetype, h5err)
            call h5tclose_f(memtype, h5err)
        end subroutine write_h5_character_scalar_attrib


        subroutine write_h5_character_vector_attrib(group, name, val)
            integer(hid_t),    intent(in)     :: group
            character(*),      intent(in)     :: name
            character(len=16), intent(in)     :: val(:)
            integer(hid_t)                    :: attr, space, filetype, memtype
            integer(hsize_t),  dimension(1:1) :: dims
            integer(size_t),   parameter      :: sdim = 16

            dims = size(val)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! type in file
            call h5tcopy_f(H5T_C_S1, filetype, h5err)
            call h5tset_size_f(filetype, sdim+1, h5err)

            ! type in run
            call h5tcopy_f(H5T_FORTRAN_S1, memtype, h5err)
            call h5tset_size_f(memtype, sdim, h5err)

            ! create the dataset
            call h5acreate_f(group, name, filetype, space, attr, h5err)

            ! write attribte
            call h5awrite_f(attr, memtype, val, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
            call h5tclose_f(filetype, h5err)
            call h5tclose_f(memtype, h5err)
        end subroutine write_h5_character_vector_attrib


        subroutine write_h5_logical_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            logical,          intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1
            integer                          :: dummy

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)

            if (val) then
                dummy = 1
            endif

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, dummy, dims, h5err)

            ! close all
            call h5aclose_f(attr, h5err)
            call h5sclose_f(space, h5err)
        end subroutine write_h5_logical_attrib

end module writer
