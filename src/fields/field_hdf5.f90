module field_hdf5
    use options, only : verbose
    use h5_utils
    use h5_writer
    use fields
    use field_diagnostics
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: hdf5_field_timer

    character(len=512) :: h5fname
    integer(hid_t)     :: h5file_id

    private :: h5file_id, h5fname

    contains

        subroutine create_h5_field_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical                  :: overwrite
            logical                  :: exists = .true.

            h5fname =  basename // '_fields.hdf5'

            ! check whether file exists
            inquire(file=h5fname, exist=exists)

            if (exists .and. overwrite) then
                call delete_h5_file(trim(h5fname))
            else if (exists) then
                print *, "File '" // trim(h5fname) // "' already exists. Exiting."
                stop
            endif

            call create_h5_file(h5fname, h5file_id)

            call write_h5_char_scalar_attrib(h5file_id, 'output_type', 'fields')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id)

            call close_h5_file(h5file_id)

        end subroutine create_h5_field_file

        subroutine write_h5_field_step(nw, t, dt)
            use options, only : output
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

            call start_timer(hdf5_field_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a18)", "write fields to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "t", t)

            call write_h5_double_scalar_step_attrib(h5file_id, nw, "dt", dt)

            call write_h5_fields(nw)

            call write_h5_field_diagnostics(h5file_id, nw)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)

            call close_h5_file(h5file_id)

            call stop_timer(hdf5_field_timer)

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

#ifndef NDEBUG
            call write_h5_dataset_2d(h5file_id, name, "symmetry volume", &
                                     sym_volg(0:nz, 0:nx-1))
#endif

            call write_h5_dataset_2d(h5file_id, name, "total buoyancy", &
                                     tbuoyg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "dry buoyancy", &
                                     dbuoyg(0:nz, 0:nx-1))

!             call write_h5_dataset_2d(h5file_id, name, "humidity", &
!                                      humg(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "vorticity", &
                                     vortg(0:nz, 0:nx-1))

!             call write_h5_dataset_2d(h5file_id, name, "liquid humidity", &
!                                      humlig(0:nz, 0:nx-1))

            call write_h5_dataset_2d(h5file_id, name, "tendency", &
                                     vtend(0:nz, 0:nx-1))

            call write_h5_int_dataset_2d(h5file_id, name, "num parcels per cell", &
                                         nparg(0:nz-1, :))

            call close_h5_group(group)
        end subroutine write_h5_fields

end module field_hdf5
