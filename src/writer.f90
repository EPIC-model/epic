!
!
! Fortran HDF5 examples at
! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html
module writer
    use hdf5
    use parcel_container, only : n_parcels, parcels
    implicit none

    private

    !character(:), allocatable :: filename

    ! h5 file handle
    integer(HID_T) :: h5file

    ! if non-zero an error occurred
    integer :: h5err = 0

    public :: h5_write_parcels, h5_create_file

    contains
        subroutine h5_create_file
            call h5fcreate_f('test.h5', H5F_ACC_TRUNC_F, h5file, h5err)

            if (h5err .ne. 0) then
                print *, "Could not create the H5 file."
                stop
            endif
        end subroutine h5_create_file

        subroutine h5_write_parcels(iter)
            integer, intent(in) :: iter ! iteration
            integer(hid_t)  :: dset ! Handles
            integer(hsize_t) :: dims(2)
            integer(hsize_t) :: space

            dims(1) = 2
            dims(2) = n_parcels
            space = n_parcels

            print *, shape(parcels%pos)

            call h5_open_file

            !
            ! Create dataspace.  Setting size to be the current size.
            !
            call h5screate_simple_f(2, dims, space, h5err)
            !
            ! Create the dataset.  We will use all default properties for this
            ! example.
            !
            call h5dcreate_f(h5file, 'step#0', H5T_NATIVE_REAL, space, dset, h5err)

            !
            ! Write the data to the dataset.
            !
            call h5dwrite_f(dset, H5T_NATIVE_REAL, parcels%pos(1:n_parcels, :), dims, h5err)

            !
            ! Close and release resources.
            !
            call h5dclose_f(dset , h5err)
            call h5sclose_f(space, h5err)

            call h5_close_file
        end subroutine h5_write_parcels


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


end module writer
