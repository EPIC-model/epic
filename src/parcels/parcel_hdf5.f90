module parcel_hdf5
    use parcel_container, only : parcels, n_parcels
    use parcel_diagnostics, only : write_h5_parcel_diagnostics
    use parcel_ellipse, only : get_angles
    use options, only : verbose, parcel, time, output, stepper
    use hdf5
    use h5_utils
    use h5_writer
    implicit none

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    private :: h5file_id, h5fname

    contains

        subroutine create_h5_parcel_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical                  :: overwrite
            logical                  :: exists = .true.

            h5fname =  basename // '_parcels.hdf5'

            ! check whether file exists
            inquire(file=h5fname, exist=exists)

            if (exists .and. overwrite) then
                call delete_h5_file(trim(h5fname))
            else if (exists) then
                print *, "File '" // trim(h5fname) // "' already exists. Exiting."
                stop
            endif

            call create_h5_file(h5fname, h5file_id)

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id)

            call close_h5_file(h5file_id)

        end subroutine create_h5_parcel_file

        subroutine write_h5_parcel_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write parcels to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "t", t)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "dt", dt)

            call write_h5_int_scalar_step_attrib(h5file_id, nw, "num parcel", n_parcels)

            call write_h5_parcel_diagnostics(h5file_id, nw)

            call write_h5_parcels(nw)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_file(h5file_id)


        end subroutine write_h5_parcel_step


        subroutine write_h5_parcels(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            double precision, allocatable :: angle(:)
            logical                       :: created

            name = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif

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
