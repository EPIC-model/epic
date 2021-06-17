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

        ! open existing group
        function open_h5_group(name) result(group)
            character(*),   intent(in) :: name
            integer(hid_t)             :: group
            logical                    :: link_exists = .false.

            call h5lexists_f(h5file, name, link_exists, h5err)

            if (link_exists) then
                call h5gopen_f(h5file, name, group, h5err)
            else
                print *, "Group '", name, "' does not exist!"
                stop
            endif
        end function open_h5_group


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

            allocate(buffer(0:dims(1)-1, 0:dims(2)-1))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)

            call h5dclose_f(dset_id, h5err)
        end subroutine read_h5_dataset_2d

        subroutine read_h5_vector_integer_attrib(group, name, buf)
            integer(hid_t), intent(in)       :: group
            character(*),   intent(in)       :: name
            integer,        intent(out)      :: buf(2)
            integer(hid_t)                   :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1) :: dims = 2

            call h5aopen_name_f(group, name, attr_id, h5err)
            call h5aget_space_f(attr_id, space_id, h5err)
            call h5aget_type_f(attr_id, type_id, h5err)
            call h5aread_f(attr_id, type_id, buf, dims, h5err)
            call h5aclose_f(attr_id, h5err)
        end subroutine read_h5_vector_integer_attrib


        subroutine read_h5_vector_double_attrib(group, name, buf)
            integer(hid_t),   intent(in)       :: group
            character(*),     intent(in)       :: name
            double precision, intent(out)      :: buf(2)
            integer(hid_t)                     :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1)   :: dims = 2

            call h5aopen_name_f(group, name, attr_id, h5err)
            call h5aget_space_f(attr_id, space_id, h5err)
            call h5aget_type_f(attr_id, type_id, h5err)
            call h5aread_f(attr_id, type_id, buf, dims, h5err)
            call h5aclose_f(attr_id, h5err)
        end subroutine read_h5_vector_double_attrib

        subroutine read_h5_box(nx, nz, extent, origin)
            integer,          intent(out)    :: nx, nz
            double precision, intent(out)    :: extent(2), origin(2)
            integer(hid_t)                   :: group
            integer                          :: ncells(2)

            group = open_h5_group("box")

            call read_h5_vector_integer_attrib(group, 'ncells', ncells)
            call read_h5_vector_double_attrib(group, 'extent', extent)
            call read_h5_vector_double_attrib(group, 'origin', origin)

            nx = ncells(1)
            nz = ncells(2)

            ! close all
            call h5gclose_f(group, h5err)
        end subroutine read_h5_box

end module reader
