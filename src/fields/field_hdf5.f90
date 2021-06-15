module field_hdf5
    use options, only : verbose
    use h5_utils
    use h5_writer
    use fields
    implicit none

    character(len=512) :: h5_field_fname
    integer(hid_t)     :: h5file_id

    private :: h5file_id

    contains

        subroutine create_h5_field_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical                  :: overwrite
            logical                  :: exists = .true.

            h5_field_fname =  basename // '_fields.hdf5'

            ! check whether file exists
            inquire(file=h5_field_fname, exist=exists)

            if (exists .and. overwrite) then
                call delete_h5_file(trim(h5_field_fname))
            else if (exists) then
                print *, "File '" // trim(h5_field_fname) // "' already exists. Exiting."
                stop
            endif

            call create_h5_file(h5_field_fname, h5file_id)
        end subroutine create_h5_field_file

        subroutine write_h5_field_step(nw, t, dt)
            use options, only : output
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a18)", "write fields to h5"
            endif
#endif

            call open_h5_file(h5_field_fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "t", t)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "dt", dt)

            call write_h5_fields(nw)

!             call write_h5_field_diagnostics(nw)

            ! update number of iterations to h5 file
!             call write_h5_num_steps(nw+1)

            call close_h5_file(h5file_id)

            ! increment counter
            nw = nw + 1

        end subroutine write_h5_field_step


        subroutine write_h5_fields(iter)
            integer, intent(in)        :: iter ! iteration
            integer(hid_t)             :: group
            character(:), allocatable  :: name
            logical                    :: created

            name = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif

            !
            ! write fields (do not write halo cells)
            !

            call write_h5_dataset_3d(h5file_id, name, "velocity", &
                                     velog(0:nz, 0:nx-1, :))

            call write_h5_dataset_3d(h5file_id, name, "velocity gradient tensor", &
                                     velgradg(0:nz, 0:nx-1, :))

            call write_h5_dataset_2d(h5file_id, name, "volume", &
                                     volg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "total buoyancy", &
                                     tbuoyg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "humidity", &
                                     humg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "vorticity", &
                                     vortg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "liquid humidity", &
                                     humlig(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "tendency", &
                                     vtend(0:nz, 0:nx-1))

            call write_h5_int_dataset_2d(h5file_id, name, "num parcels per cell", &
                                         nparg(0:nz-1, :))

            call close_h5_group(group)
        end subroutine write_h5_fields

end module field_hdf5
