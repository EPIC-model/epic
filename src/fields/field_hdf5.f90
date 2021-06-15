module field_hdf5
    use options, only : verbose
    use h5_utils
    use writer, only : write_h5_double_scalar_step_attrib,  &
                       write_h5_dataset_1d,                 &
                       write_h5_dataset_2d,                 &
                       write_h5_int_dataset_2d,             &
                       write_h5_dataset_3d,                 &
                       get_step_group_name
    use fields
    implicit none

    private
        character(len=512) :: h5fname
        integer(hid_t)     :: h5file_id

    public :: write_h5_field_step


    contains

        subroutine write_h5_field_step(nw, t, dt)
            use options, only : output
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            logical                         :: exists = .false.

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a30)", "write fields to h5"
            endif
#endif

            ! check whether file exists
            inquire(file=h5fname, exist=exists)
            if (exists .eqv. .false.) then
                call create_h5_file(h5fname, h5file_id)
            endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

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

            name = trim(get_step_group_name(iter))

            call open_or_create_h5_group(h5file_id, name, group)

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
