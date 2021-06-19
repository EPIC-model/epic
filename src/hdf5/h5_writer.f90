! =============================================================================
! This module is used to write data into a HDF5 file. Examples of HDF5 can be
! found at https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! The API specifications are found in
! https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5Front.html
! =============================================================================
module h5_writer
    use options, only : allow_larger_anisotropy,    &
                        output, verbose,            &
                        field_file, field_tol,      &
                        parcel, time
    use parameters, only : nx, nz, lower, extent
    use hdf5
    use h5_utils
    implicit none

    contains

        subroutine write_h5_dataset_1d(h5file_id, group, name, data)
            integer(hid_t),   intent(in)     :: h5file_id
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:1) :: dims

            if (size(data) == 0) then
                print *, "Error in 'write_h5_dataset_1d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(1, dims, dataspace, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5dcreate_f(h5file_id, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)
            call check_h5_error("Failed to create dataset.")

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)
            call check_h5_error("Failed to write dataset.")

            ! close all
            call h5dclose_f(dset , h5err)
            call check_h5_error("Failed to close dataset.")

            call h5sclose_f(dataspace, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_dataset_1d

        subroutine write_h5_dataset_2d(h5file_id, group, name, data)
            integer(hid_t),   intent(in)     :: h5file_id
            ! 12 March 2021
            ! https://stackoverflow.com/questions/48816383/passing-character-strings-of-different-lengths-to-functions-in-fortran
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:, :)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:2) :: dims

            if (size(data) == 0) then
                print *, "Error in 'write_h5_dataset_2d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(2, dims, dataspace, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5dcreate_f(h5file_id, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)
            call check_h5_error("Failed to create dataset.")

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)
            call check_h5_error("Failed to write dataset.")

            ! close all
            call h5dclose_f(dset , h5err)
            call check_h5_error("Failed to close dataset.")

            call h5sclose_f(dataspace, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_dataset_2d

        subroutine write_h5_int_dataset_2d(h5file_id, group, name, data)
            integer(hid_t),   intent(in)     :: h5file_id
            ! 12 March 2021
            ! https://stackoverflow.com/questions/48816383/passing-character-strings-of-different-lengths-to-functions-in-fortran
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: data(:, :)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:2) :: dims

            if (size(data) == 0) then
                print *, "Error in 'write_h5_int_dataset_2d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(2, dims, dataspace, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5dcreate_f(h5file_id, group // "/" // name,            &
                             H5T_NATIVE_INTEGER, dataspace, dset, h5err)
            call check_h5_error("Failed to create dataset.")

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_INTEGER, data, dims, h5err)
            call check_h5_error("Failed to write dataset.")

            ! close all
            call h5dclose_f(dset , h5err)
            call check_h5_error("Failed to close dataset.")

            call h5sclose_f(dataspace, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_int_dataset_2d

        subroutine write_h5_dataset_3d(h5file_id, group, name, data)
            integer(hid_t),   intent(in)     :: h5file_id
            character(*),     intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: data(:, :, :)
            integer(hid_t)                   :: dset, dataspace
            integer(hsize_t), dimension(1:3) :: dims

            if (size(data) == 0) then
                print *, "Error in 'write_h5_dataset_3d': ", &
                         "No memory for '", name, "' allocated!"
                stop
            endif

            dims = shape(data)

            ! create space for data
            call h5screate_simple_f(3, dims, dataspace, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5dcreate_f(h5file_id, group // "/" // name,              &
                             H5T_NATIVE_DOUBLE, dataspace, dset, h5err)
            call check_h5_error("Failed to create dataset.")

            ! write dataset
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, data, dims, h5err)
            call check_h5_error("Failed to write dataset.")

            ! close all
            call h5dclose_f(dset , h5err)
            call check_h5_error("Failed to close dataset.")

            call h5sclose_f(dataspace, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_dataset_3d


        subroutine write_h5_double_scalar_step_attrib(h5file_id, iter, name, val)
            integer(hid_t),   intent(in)     :: h5file_id
            integer,          intent(in)     :: iter ! iteration
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: val
            integer(hid_t)                   :: group
            character(:), allocatable        :: grn
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1
            logical                          :: created

            grn = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, grn, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, grn, group)
            endif

            call write_h5_double_scalar_attrib(group, name, val)

            call close_h5_group(group)
        end subroutine write_h5_double_scalar_step_attrib


        subroutine write_h5_int_scalar_step_attrib(h5file_id, iter, name, val)
            integer(hid_t),   intent(in)     :: h5file_id
            integer,          intent(in)     :: iter ! iteration
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: val
            integer(hid_t)                   :: group
            character(:), allocatable        :: grn
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1
            logical                          :: created

            grn = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, grn, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, grn, group)
            endif

            call write_h5_int_scalar_attrib(group, name, val)

            call close_h5_group(group)
        end subroutine write_h5_int_scalar_step_attrib


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
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_DOUBLE, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            call h5awrite_f(attr, H5T_NATIVE_DOUBLE, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute space.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_double_vector_attrib

        subroutine write_h5_double_scalar_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            double precision, intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_DOUBLE, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            call h5awrite_f(attr, H5T_NATIVE_DOUBLE, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute space.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_double_scalar_attrib

        subroutine write_h5_int_scalar_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute space.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_int_scalar_attrib


        subroutine write_h5_int_vector_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            integer,          intent(in)     :: val(:)
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims

            dims = size(val)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute space.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_int_vector_attrib


        subroutine write_h5_char_scalar_attrib(group, name, val)
            integer(hid_t),    intent(in)     :: group
            character(*),      intent(in)     :: name
            character(len=*), intent(in)     :: val
            integer(hid_t)                    :: attr, space, filetype, memtype
            integer(hsize_t),  dimension(1:1) :: dims = 1
            integer(size_t)                   :: sdim

            sdim = len(val)

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! type in file
            call h5tcopy_f(H5T_C_S1, filetype, h5err)
            call h5tset_size_f(filetype, sdim+1, h5err)

            ! type in run
            call h5tcopy_f(H5T_FORTRAN_S1, memtype, h5err)
            call check_h5_error("Failed to copy memory type.")

            call h5tset_size_f(memtype, sdim, h5err)
            call check_h5_error("Failed to create memory space.")

            ! create the dataset
            call h5acreate_f(group, name, filetype, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            ! write attribte
            call h5awrite_f(attr, memtype, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute space.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")

            call h5tclose_f(filetype, h5err)
            call check_h5_error("Failed to close file type space.")

            call h5tclose_f(memtype, h5err)
            call check_h5_error("Failed to close memory space.")
        end subroutine write_h5_char_scalar_attrib


        subroutine write_h5_char_vector_attrib(group, name, val)
            integer(hid_t),    intent(in)     :: group
            character(*),      intent(in)     :: name
            character(len=*), intent(in)      :: val(:)
            integer(hid_t)                    :: attr, space, filetype, memtype
            integer(hsize_t),  dimension(1:1) :: dims
            integer(size_t)                   :: sdim

            dims = size(val)
            sdim = len(val(1))

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! type in file
            call h5tcopy_f(H5T_C_S1, filetype, h5err)
            call check_h5_error("Failed to copy file type.")

            call h5tset_size_f(filetype, sdim+1, h5err)
            call check_h5_error("Failed to set size.")

            ! type in run
            call h5tcopy_f(H5T_FORTRAN_S1, memtype, h5err)
            call check_h5_error("Failed to copy memory sace.")

            call h5tset_size_f(memtype, sdim, h5err)
            call check_h5_error("Failed to set memory size.")

            ! create the dataset
            call h5acreate_f(group, name, filetype, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            ! write attribte
            call h5awrite_f(attr, memtype, val, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")

            call h5tclose_f(filetype, h5err)
            call check_h5_error("Failed to file type space.")

            call h5tclose_f(memtype, h5err)
            call check_h5_error("Failed to close memory space.")
        end subroutine write_h5_char_vector_attrib


        subroutine write_h5_logical_attrib(group, name, val)
            integer(hid_t),   intent(in)     :: group
            character(*),     intent(in)     :: name
            logical,          intent(in)     :: val
            integer(hid_t)                   :: attr, space
            integer(hsize_t), dimension(1:1) :: dims = 1
            integer                          :: dummy

            ! create space for data
            call h5screate_simple_f(1, dims, space, h5err)
            call check_h5_error("Failed to create data space.")

            ! create the dataset
            call h5acreate_f(group, name,                           &
                             H5T_NATIVE_INTEGER, space, attr, h5err)
            call check_h5_error("Failed to create attribute.")

            dummy = 0
            if (val) then
                dummy = 1
            endif

            call h5awrite_f(attr, H5T_NATIVE_INTEGER, dummy, dims, h5err)
            call check_h5_error("Failed to write attribute.")

            ! close all
            call h5aclose_f(attr, h5err)
            call check_h5_error("Failed to close attribute.")

            call h5sclose_f(space, h5err)
            call check_h5_error("Failed to close data space.")
        end subroutine write_h5_logical_attrib

        subroutine write_h5_num_steps(h5file_id, nw)
            integer(hid_t),   intent(in) :: h5file_id
            integer,          intent(in) :: nw
            integer(hid_t)               :: group
            logical                      :: attr_exists

            call open_h5_group(h5file_id, "/", group)

            ! in principle not necessary but we still check
            call h5aexists_f(group, "nsteps", attr_exists, h5err)
            call check_h5_error("Failed to check if attribute exists.")

            if (attr_exists) then
                call h5adelete_f(group, "nsteps", h5err)
                call check_h5_error("Failed to delete attribute.")
            endif

            call write_h5_int_scalar_attrib(group, "nsteps", nw)

            call close_h5_group(group)

        end subroutine write_h5_num_steps

        subroutine write_h5_timestamp(h5file_id)
            integer(hid_t),   intent(in) :: h5file_id
            character(len=8)             :: cdate, tmp2
            character(len=10)            :: ctime, tmp1
            character(len=5)             :: czone
            character(len=9)             :: tmp3

            ! 15 June 2021
            ! https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005fTIME.html
            call date_and_time(date=cdate, time=ctime, zone=czone)
            ! 15 June 2021
            ! https://stackoverflow.com/questions/13755762/access-character-at-specific-index-in-a-string-in-fortran
            tmp1 = cdate(1:4) // '/' // cdate(5:6) // '/' // cdate(7:8)
            tmp2 = ctime(1:2) // ':' // ctime(3:4) // ':' // ctime(5:6)
            tmp3 = 'UTC' // czone(1:3) // ':' // czone(4:5)
            call write_h5_char_scalar_attrib(h5file_id, "creation date", tmp1)
            call write_h5_char_scalar_attrib(h5file_id, "creation time", tmp2)
            call write_h5_char_scalar_attrib(h5file_id, "creation zone", tmp3)
        end subroutine write_h5_timestamp


        subroutine write_h5_options(h5file_id)
            integer(hid_t),   intent(in) :: h5file_id
            integer(hid_t)               :: gopts, group

            call create_h5_group(h5file_id, "options", gopts)

#ifdef ENABLE_VERBOSE
            call write_h5_logical_attrib(gopts, "verbose", verbose)
#endif
            call write_h5_char_scalar_attrib(gopts, "field_file", field_file)
            call write_h5_double_scalar_attrib(gopts, "field_tol", field_tol)

            call create_h5_group(gopts, "parcel", group)
                call write_h5_int_scalar_attrib(group, "n_per_cell", parcel%n_per_cell)
                call write_h5_logical_attrib(group, "is_elliptic", parcel%is_elliptic)
                call write_h5_double_scalar_attrib(group, "lambda", parcel%lambda)
                call write_h5_double_scalar_attrib(group, "prefactor", parcel%prefactor)
                call write_h5_int_scalar_attrib(group, "split_freq", parcel%split_freq)
                call write_h5_char_scalar_attrib(group, "merge_type", trim(parcel%merge_type))
                call write_h5_int_scalar_attrib(group, "merge_freq", parcel%merge_freq)
                call write_h5_double_scalar_attrib(group, "vfraction", parcel%vfraction)
                call write_h5_int_scalar_attrib(group, "correction_freq", parcel%correction_freq)
                call write_h5_int_scalar_attrib(group, "correction_iters", parcel%correction_iters)
                call write_h5_logical_attrib(group, "apply_laplace", parcel%apply_laplace)
                call write_h5_logical_attrib(group, "apply_gradient", parcel%apply_gradient)
                call write_h5_double_scalar_attrib(group, "gradient_pref", parcel%gradient_pref)
                call write_h5_double_scalar_attrib(group, "vmaxfraction", parcel%vmaxfraction)
            call close_h5_group(group)

            call create_h5_group(gopts, "output", group)
                call write_h5_int_scalar_attrib(group, "h5_parcel_freq", output%h5_parcel_freq)
                call write_h5_int_scalar_attrib(group, "h5_field_freq", output%h5_field_freq)
                call write_h5_logical_attrib(group, "h5_write_fields", output%h5_write_fields)
                call write_h5_logical_attrib(group, "h5_overwrite", output%h5_overwrite)
                call write_h5_logical_attrib(group, "h5_write_parcels", output%h5_write_parcels)
                call write_h5_char_scalar_attrib(group, "h5_basename", trim(output%h5_basename))
            call close_h5_group(group)

            call create_h5_group(gopts, "time", group)
                call write_h5_double_scalar_attrib(group, "limit", time%limit)
                call write_h5_logical_attrib(group, "is_adaptive", time%is_adaptive)
                call write_h5_double_scalar_attrib(group, "alpha", time%alpha)
                call write_h5_double_scalar_attrib(group, "dt", time%dt)
                call write_h5_double_scalar_attrib(group, "dt_max", time%dt_max)
            call close_h5_group(group)

            call close_h5_group(gopts)
        end subroutine write_h5_options

        subroutine write_h5_box(h5file_id)
            integer(hid_t), intent(in) :: h5file_id
            integer(hid_t)             :: group

            call create_h5_group(h5file_id, "box", group)
            call write_h5_logical_attrib(group, "allow_larger_anisotropy", &
                                         allow_larger_anisotropy)
            call write_h5_int_vector_attrib(group, "ncells", (/nx, nz/))
            call write_h5_double_vector_attrib(group, "extent", extent)
            call write_h5_double_vector_attrib(group, "origin", lower)
            call close_h5_group(group)
        end subroutine write_h5_box

end module h5_writer
