! =============================================================================
!                               HDF5 reader
!   This module reads in datatsets. For further information about the
!   HDF5 API please visit
!   https://support.hdfgroup.org/HDF5/doc1.6/RM_H5Front.html
! =============================================================================
module h5_reader
    use hdf5
    use h5_utils
    implicit none

    private :: read_h5_dataset_1d, &
               read_h5_dataset_2d, &
               read_h5_dataset_3d, &
               read_h5_dataset_4d

    interface read_h5_dataset
        module procedure :: read_h5_dataset_1d
        module procedure :: read_h5_dataset_2d
        module procedure :: read_h5_dataset_3d
        module procedure :: read_h5_dataset_4d
    end interface read_h5_dataset

    contains

        ! @returns -1 if H5 file does not contain the attribute 'num parcel'
        subroutine get_num_parcels(h5file_id, n_parcels)
            integer(hid_t),   intent(in)     :: h5file_id
            integer,          intent(out)    :: n_parcels
            integer(hid_t)                   :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1) :: dims = 1

            if (.not. has_attribute(h5file_id, 'num parcel')) then
                n_parcels = -1
                return
            endif

            call h5aopen_name_f(h5file_id, 'num parcel', attr_id, h5err)
            call check_h5_error("Failed to open attribute.")
            call h5aget_space_f(attr_id, space_id, h5err)
            call check_h5_error("Failed to get attribute space.")
            call h5aget_type_f(attr_id, type_id, h5err)
            call check_h5_error("Failed to retrieve attribute type.")
            call h5aread_f(attr_id, type_id, n_parcels, dims, h5err)
            call check_h5_error("Failed to read attribute.")
            call h5aclose_f(attr_id, h5err)
            call check_h5_error("Failed to close attribute.")
        end subroutine get_num_parcels

        ! @returns -1 if H5 file does not contain the attribute 'nsteps'
        subroutine get_num_steps(h5file_id, n_steps)
            integer(hid_t),   intent(in)     :: h5file_id
            integer,          intent(out)    :: n_steps
            integer(hid_t)                   :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1) :: dims = 1

            if (.not. has_attribute(h5file_id, 'nsteps')) then
                n_steps = -1
                return
            endif

            call h5aopen_name_f(h5file_id, 'nsteps', attr_id, h5err)
            call check_h5_error("Failed to open attribute.")
            call h5aget_space_f(attr_id, space_id, h5err)
            call check_h5_error("Failed to get attribute space.")
            call h5aget_type_f(attr_id, type_id, h5err)
            call check_h5_error("Failed to retrieve attribute type.")
            call h5aread_f(attr_id, type_id, n_steps, dims, h5err)
            call check_h5_error("Failed to read attribute.")
            call h5aclose_f(attr_id, h5err)
            call check_h5_error("Failed to close attribute.")
        end subroutine get_num_steps

        ! 11 Jan 2022
        ! https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/FORTRAN/H5T/h5ex_t_stringCatt_F03.f90
        subroutine get_file_type(h5file_id, file_type)
            use iso_c_binding
            integer(hid_t),            intent(in)  :: h5file_id
            character(:), allocatable, intent(out) :: file_type
            character(len=512), target             :: buf
            integer(hid_t)                         :: attr_id, space_id, type_id, memtype
            integer(hsize_t)                       :: sdim
            integer(hsize_t)                       :: length
            type(c_ptr)                            :: f_ptr

            if (.not. has_attribute(h5file_id, 'output_type')) then
                print *, 'Not a proper EPIC HDF5 file.'
                stop
            endif

            call h5aopen_f(h5file_id, 'output_type', attr_id, h5err)
            call check_h5_error("Failed to open attribute.")

            call h5aget_type_f(attr_id, type_id, h5err)
            call check_h5_error("Failed to retrieve attribute type.")

            call h5tget_size_f(type_id, length, h5err)
            call check_h5_error("Failed to get size.")

            sdim = length

            call h5aget_space_f(attr_id, space_id, h5err)

            call h5tcopy_f(H5T_FORTRAN_S1, memtype, h5err)
            call h5tset_size_f(memtype, sdim, h5err)

            f_ptr = c_loc(buf)
            call h5aread_f(attr_id, memtype, f_ptr, h5err)
            call check_h5_error("Failed to read attribute.")
            call h5aclose_f(attr_id, h5err)
            call check_h5_error("Failed to close attribute.")

            file_type = buf(1:length)
        end subroutine get_file_type

        function has_attribute(h5file_id, name) result(link_exists)
            integer(hid_t), intent(in) :: h5file_id
            character(*),   intent(in) :: name
            logical                    :: link_exists
            link_exists = .false.
            call h5aexists_f(h5file_id, name, link_exists, h5err)
            call check_h5_error("Failed to check if link exists.")
        end function has_attribute

        function has_dataset(h5file_id, name) result(link_exists)
            integer(hid_t), intent(in) :: h5file_id
            character(*),   intent(in) :: name
            logical                    :: link_exists
            link_exists = .false.
            call h5lexists_f(h5file_id, name, link_exists, h5err)
            call check_h5_error("Failed to check if link exists.")
        end function has_dataset

        subroutine read_h5_dataset_1d(h5file_id, name, buffer)
            integer(hid_t),                intent(in)  :: h5file_id
            character(*),                  intent(in)  :: name
            double precision, allocatable, intent(out) :: buffer(:)
            integer(hid_t)                             :: dataspace_id   ! Dataspace identifier
            integer(hid_t)                             :: dset_id
            integer                                    :: rank
            integer(hsize_t)                           :: dims(1), maxdims(1)

            call h5dopen_f(h5file_id, name, dset_id, h5err)
            call check_h5_error("Failed to open dataset.")

            call h5dget_space_f(dset_id, dataspace_id, h5err)
            call check_h5_error("Failed to get data space.")

            call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)
            call check_h5_error("Failed to get data extent.")

            if (rank .ne. 1) then
                print *, "Dataset not 1-dimensional."
                stop
            endif

            ! h5err: Dataspace rank on success and -1 on failure
            call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, h5err)
#ifndef NDEBUG
            if (h5err >= 0) then
                h5err = 0
            endif
#endif
            call check_h5_error("Failed to get dimensions.")

            if (.not. sum(dims -  maxdims) == 0) then
                print *, "Dimensions do not agree."
                stop
            endif

            allocate(buffer(0:dims(1)-1))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)
            call check_h5_error("Failed to read dataset.")

            call h5dclose_f(dset_id, h5err)
            call check_h5_error("Failed to close dataset.")
        end subroutine read_h5_dataset_1d

        subroutine read_h5_dataset_2d(h5file_id, name, buffer)
            integer(hid_t),                intent(in)  :: h5file_id
            character(*),                  intent(in)  :: name
            double precision, allocatable, intent(out) :: buffer(:, :)
            integer(hid_t)                             :: dataspace_id   ! Dataspace identifier
            integer(hid_t)                             :: dset_id
            integer                                    :: rank
            integer(hsize_t)                           :: dims(2), maxdims(2)

            call h5dopen_f(h5file_id, name, dset_id, h5err)
            call check_h5_error("Failed to open dataset.")

            call h5dget_space_f(dset_id, dataspace_id, h5err)
            call check_h5_error("Failed to get data space.")

            call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)
            call check_h5_error("Failed to get data extent.")

            if (rank .ne. 2) then
                print *, "Dataset not 2-dimensional."
                stop
            endif

            ! h5err: Dataspace rank on success and -1 on failure
            call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, h5err)
#ifndef NDEBUG
            if (h5err >= 0) then
                h5err = 0
            endif
#endif
            call check_h5_error("Failed to get dimensions.")

            if (.not. sum(dims -  maxdims) == 0) then
                print *, "Dimensions do not agree."
                stop
            endif

            allocate(buffer(0:dims(1)-1, 0:dims(2)-1))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)
            call check_h5_error("Failed to read dataset.")

            call h5dclose_f(dset_id, h5err)
            call check_h5_error("Failed to close dataset.")
        end subroutine read_h5_dataset_2d

        subroutine read_h5_dataset_3d(h5file_id, name, buffer)
            integer(hid_t),                intent(in)  :: h5file_id
            character(*),                  intent(in)  :: name
            double precision, allocatable, intent(out) :: buffer(:, :, :)
            integer(hid_t)                             :: dataspace_id   ! Dataspace identifier
            integer(hid_t)                             :: dset_id
            integer                                    :: rank
            integer(hsize_t)                           :: dims(3), maxdims(3)

            call h5dopen_f(h5file_id, name, dset_id, h5err)
            call check_h5_error("Failed to open dataset.")

            call h5dget_space_f(dset_id, dataspace_id, h5err)
            call check_h5_error("Failed to get data space.")

            call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)
            call check_h5_error("Failed to get data extent.")

            if (rank .ne. 3) then
                print *, "Dataset not 3-dimensional."
                stop
            endif

            ! h5err: Dataspace rank on success and -1 on failure
            call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, h5err)
#ifndef NDEBUG
            if (h5err >= 0) then
                h5err = 0
            endif
#endif
            call check_h5_error("Failed to get dimensions.")

            if (.not. sum(dims -  maxdims) == 0) then
                print *, "Dimensions do not agree."
                stop
            endif

            allocate(buffer(0:dims(1)-1, 0:dims(2)-1, 0:dims(3)-1))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)
            call check_h5_error("Failed to read dataset.")

            call h5dclose_f(dset_id, h5err)
            call check_h5_error("Failed to close dataset.")
        end subroutine read_h5_dataset_3d

        subroutine read_h5_dataset_4d(h5file_id, name, buffer)
            integer(hid_t),                intent(in)  :: h5file_id
            character(*),                  intent(in)  :: name
            double precision, allocatable, intent(out) :: buffer(:, :, :, :)
            integer(hid_t)                             :: dataspace_id   ! Dataspace identifier
            integer(hid_t)                             :: dset_id
            integer                                    :: rank
            integer(hsize_t)                           :: dims(4), maxdims(4)

            call h5dopen_f(h5file_id, name, dset_id, h5err)
            call check_h5_error("Failed to open dataset.")

            call h5dget_space_f(dset_id, dataspace_id, h5err)
            call check_h5_error("Failed to get data space.")

            call h5sget_simple_extent_ndims_f(dataspace_id, rank, h5err)
            call check_h5_error("Failed to get data extent.")

            if (rank .ne. 4) then
                print *, "Dataset not 4-dimensional."
                stop
            endif

            ! h5err: Dataspace rank on success and -1 on failure
            call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, h5err)
#ifndef NDEBUG
            if (h5err >= 0) then
                h5err = 0
            endif
#endif
            call check_h5_error("Failed to get dimensions.")

            if (.not. sum(dims -  maxdims) == 0) then
                print *, "Dimensions do not agree."
                stop
            endif

            allocate(buffer(0:dims(1)-1, 0:dims(2)-1, 0:dims(3)-1, 1:dims(4)))

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, h5err)
            call check_h5_error("Failed to read dataset.")

            call h5dclose_f(dset_id, h5err)
            call check_h5_error("Failed to close dataset.")
        end subroutine read_h5_dataset_4d

        subroutine read_h5_vector_int_attrib(group, name, buf)
            integer(hid_t), intent(in)       :: group
            character(*),   intent(in)       :: name
            integer,        intent(out)      :: buf(2)
            integer(hid_t)                   :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1) :: dims = 2

            call h5aopen_name_f(group, name, attr_id, h5err)
            call check_h5_error("Failed to open attribute.")
            call h5aget_space_f(attr_id, space_id, h5err)
            call check_h5_error("Failed to get attribute space.")
            call h5aget_type_f(attr_id, type_id, h5err)
            call check_h5_error("Failed to retrieve attribute type.")
            call h5aread_f(attr_id, type_id, buf, dims, h5err)
            call check_h5_error("Failed to read attribute.")
            call h5aclose_f(attr_id, h5err)
            call check_h5_error("Failed to close attribute.")
        end subroutine read_h5_vector_int_attrib


        subroutine read_h5_vector_double_attrib(group, name, buf)
            integer(hid_t),   intent(in)       :: group
            character(*),     intent(in)       :: name
            double precision, intent(out)      :: buf(2)
            integer(hid_t)                     :: attr_id, space_id, type_id
            integer(hsize_t), dimension(1:1)   :: dims = 2

            call h5aopen_name_f(group, name, attr_id, h5err)
            call check_h5_error("Failed to open attribute.")
            call h5aget_space_f(attr_id, space_id, h5err)
            call check_h5_error("Failed to get attribute space.")
            call h5aget_type_f(attr_id, type_id, h5err)
            call check_h5_error("Failed to retrieve attribute type.")
            call h5aread_f(attr_id, type_id, buf, dims, h5err)
            call check_h5_error("Failed to read attribute.")
            call h5aclose_f(attr_id, h5err)
            call check_h5_error("Failed to close attribute.")
        end subroutine read_h5_vector_double_attrib

        subroutine read_h5_box(h5file_id, ncells, extent, origin)
            integer(hid_t),   intent(in)     :: h5file_id
            integer,          intent(out)    :: ncells(:)
            double precision, intent(out)    :: extent(:), origin(:)
            integer(hid_t)                   :: group

            if ((size(ncells) > 3) .or. (size(extent) > 3) .or. (size(extent) > 3)) then
                print *, "Cannot read more than 3 dimensions!"
                stop
            endif

            call open_h5_group(h5file_id, "box", group)

            call read_h5_vector_int_attrib(group, 'ncells', ncells)
            call read_h5_vector_double_attrib(group, 'extent', extent)
            call read_h5_vector_double_attrib(group, 'origin', origin)

            ! close all
            call close_h5_group(group)
        end subroutine read_h5_box

end module h5_reader
