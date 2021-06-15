module parcel_hdf5
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics, only : write_h5_parcel_diagnostics
    use parcel_ellipse, only : get_angles
    use options, only : verbose
    use hdf5
    use h5_utils
    use writer, only : write_h5_double_scalar_step_attrib,  &
                       write_h5_int_scalar_step_attrib, &
                       write_h5_int_scalar_attrib,      &
                       write_h5_dataset_1d, &
                       write_h5_dataset_2d, &
                       get_step_group_name
    implicit none

    private
        character(len=512) :: h5fname

        ! h5 file handle
        integer(hid_t) :: h5file_id


    public :: write_h5_parcel_step


    contains

        subroutine write_h5_parcel_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            logical                         :: exists = .false.

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a30)", "write parcels to h5"
            endif
#endif
            ! check whether file exists
            inquire(file=h5fname, exist=exists)
            if (exists .eqv. .false.) then
                call create_h5_file(h5fname, h5file_id)
            endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_int_scalar_step_attrib(h5file_id, nw, "num parcel", n_parcels)

            call write_h5_parcels(nw)

            call write_h5_parcel_diagnostics(h5file_id, nw)

            ! update number of iterations to h5 file
!             call write_h5_num_steps(nw+1)

            call close_h5_file(h5file_id)

            ! increment counter
            nw = nw + 1

        end subroutine write_h5_parcel_step


        subroutine write_h5_parcels(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            double precision, allocatable :: angle(:)

            name = trim(get_step_group_name(iter))

            call open_or_create_h5_group(h5file_id, name, group)

            !
            ! write parcel data
            !

            call write_h5_dataset_2d(h5file_id, name, "position", &
                                     parcels%position(1:n_parcels, :))

            call write_h5_dataset_2d(h5file_id, name, "velocity", &
                                     parcels%velocity(1:n_parcels, :))

            call write_h5_dataset_1d(h5file_id, name, "vorticity", &
                                     parcels%vorticity(1:n_parcels))

            if (allocated(parcels%stretch)) then
                call write_h5_dataset_1d(h5file_id, name, "stretch", &
                                         parcels%stretch(1:n_parcels))
            endif

            call write_h5_dataset_1d(h5file_id, name, "volume", &
                                     parcels%volume(1:n_parcels))

            call write_h5_dataset_1d(h5file_id, name, "buoyancy", &
                                     parcels%buoyancy(1:n_parcels))

            call write_h5_dataset_1d(h5file_id, name, "humidity", &
                                     parcels%humidity(1:n_parcels))

            if (allocated(parcels%B)) then
                call write_h5_dataset_2d(h5file_id, name, "B", &
                                         parcels%B(1:n_parcels, :))

                angle = get_angles(parcels%B, parcels%volume, n_parcels)
                call write_h5_dataset_1d(h5file_id, name, "orientation", angle)
            endif

            call close_h5_group(group)

        end subroutine write_h5_parcels

end module parcel_hdf5
