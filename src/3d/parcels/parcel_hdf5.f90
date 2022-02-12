module parcel_hdf5
    use constants, only : max_num_parcels
    use parcel_container, only : parcels, n_parcels
    use options, only : verbose, write_h5_options
    use parameters, only : nx, ny, nz, extent, lower
    use hdf5
    use h5_utils
    use h5_writer
    use h5_reader
    implicit none

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    private :: h5file_id, h5fname

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_h5_parcel_file(basename, overwrite, l_restart, step)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            integer,      intent(out) :: step
            logical                   :: l_exist

            h5fname =  basename // '_parcels.hdf5'

            call exist_h5_file(h5fname, l_exist)

            if (l_restart .and. l_exist) then
                call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)
                call get_num_steps(h5file_id, step)
                call close_h5_file(h5file_id)
                return
            endif

            step = 0

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_scalar_attrib(h5file_id, 'output_type', 'parcels')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id, lower, extent, (/nx, ny, nz/))

            call close_h5_file(h5file_id)

        end subroutine create_h5_parcel_file

        ! Write diagnostics for the current time step in the parcel file.
        ! @param[inout] nw counts the number of writes
        ! @param[in] t is the time
        ! @param[in] dt is the time step
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

            call write_h5_scalar_step_attrib(h5file_id, nw, "t", t)

            call write_h5_scalar_step_attrib(h5file_id, nw, "dt", dt)

            call write_h5_scalar_step_attrib(h5file_id, nw, "num parcel", n_parcels)

            call write_h5_parcels(nw)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_file(h5file_id)

        end subroutine write_h5_parcel_step


        ! Write parcel datasets (called from write_h5_parcel_step).
        ! @param[in] iter is the number of the write
        subroutine write_h5_parcels(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            logical                       :: created

            name = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif

            !
            ! write parcel data
            !

            call write_h5_dataset(h5file_id, name, "position", &
                                  parcels%position(:, 1:n_parcels))

            call write_h5_dataset(h5file_id, name, "B", &
                                  parcels%B(:, 1:n_parcels))

            call write_h5_dataset(h5file_id, name, "volume", &
                                  parcels%volume(1:n_parcels))

            call write_h5_dataset(h5file_id, name, "vorticity", &
                                  parcels%vorticity(:, 1:n_parcels))

            call write_h5_dataset(h5file_id, name, "buoyancy", &
                                  parcels%buoyancy(1:n_parcels))
#ifndef ENABLE_DRY_MODE
            call write_h5_dataset(h5file_id, name, "humidity", &
                                  parcels%humidity(1:n_parcels))
#endif
            call close_h5_group(group)

        end subroutine write_h5_parcels


        ! This subroutine reads parcels from an EPIC output file
        subroutine read_h5_parcels(h5fname, step)
            character(*),     intent(in)  :: h5fname
            integer,          intent(in)  :: step
            integer(hid_t)                :: h5handle, group
            integer                       :: ncells(3)
            character(:), allocatable     :: grn
            double precision, allocatable :: buffer_1d(:),   &
                                             buffer_2d(:, :)
            logical                       :: l_valid = .false.


            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)

            ! read domain dimensions
            call read_h5_box(h5handle, ncells, extent, lower)
            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)

            ! update global parameters
            call update_parameters

            grn = trim(get_step_group_name(step))
            call open_h5_group(h5handle, grn, group)
            call get_num_parcels(group, n_parcels)

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                stop
            endif

            ! Be aware that the starting index of buffer_1d and buffer_2d
            ! is 0; hence, the range is 0:n_parcels-1 in contrast to the
            ! parcel container where it is 1:n_parcels.

            if (has_dataset(group, 'B')) then
                call read_h5_dataset(group, 'B', buffer_2d)
                parcels%B(:, 1:n_parcels) = buffer_2d
                deallocate(buffer_2d)
            else
                print *, "The parcel shape must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'position')) then
                call read_h5_dataset(group, 'position', buffer_2d)
                parcels%position(:, 1:n_parcels) = buffer_2d
                deallocate(buffer_2d)
            else
                print *, "The parcel position must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'volume')) then
                call read_h5_dataset(group, 'volume', buffer_1d)
                parcels%volume(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            else
                print *, "The parcel volume must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'vorticity')) then
                l_valid = .true.
                call read_h5_dataset(group, 'vorticity', buffer_2d)
                parcels%vorticity(:, 1:n_parcels) = buffer_2d
                deallocate(buffer_2d)
            endif

            if (has_dataset(group, 'buoyancy')) then
                l_valid = .true.
                call read_h5_dataset(group, 'buoyancy', buffer_1d)
                parcels%buoyancy(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(group, 'humidity')) then
                l_valid = .true.
                call read_h5_dataset(group, 'humidity', buffer_1d)
                parcels%buoyancy(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            endif
#endif

            if (.not. l_valid) then
                print *, "Either the parcel buoyancy or vorticity must be present! Exiting."
                stop
            endif

            call close_h5_group(group)
            call close_h5_file(h5handle)

        end subroutine read_h5_parcels

end module parcel_hdf5
