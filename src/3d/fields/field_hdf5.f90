module field_hdf5
    use options, only : verbose, write_h5_options
    use h5_utils
    use h5_writer
    use h5_reader, only : get_num_steps
    use fields
    implicit none

    integer :: n_writes

    character(len=512) :: h5fname
    integer(hid_t)     :: h5file_id

    private :: h5file_id, h5fname, n_writes

    contains

        ! Create the field file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_h5_field_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist

            h5fname =  basename // '_fields.hdf5'

            call exist_h5_file(h5fname, l_exist)

            if (l_restart .and. l_exist) then
                call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)
                call get_num_steps(h5file_id, n_writes)
                call close_h5_file(h5file_id)
                return
            endif

            n_writes = 0

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_scalar_attrib(h5file_id, 'output_type', 'fields')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id, lower, extent, (/nx, ny, nz/))

            call close_h5_file(h5file_id)

        end subroutine create_h5_field_file

        ! Write a step in the field file.
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_h5_field_step(t, dt)
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a18)", "write fields to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            call write_h5_scalar_step_attrib(h5file_id, n_writes, "t", t)

            call write_h5_scalar_step_attrib(h5file_id, n_writes, "dt", dt)

            call write_h5_fields(n_writes)

            ! increment counter
            n_writes = n_writes + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, n_writes)

            call close_h5_file(h5file_id)

        end subroutine write_h5_field_step


        ! Write field datasets (called from write_h5_field_step).
        ! @param[in] iter is the number of the write
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
            call write_h5_dataset(h5file_id, name, "velocity", &
                                  velog(0:nz, 0:ny-1, 0:nx-1, :))

            call write_h5_dataset(h5file_id, name, "vorticity", &
                                  vortg(0:nz, 0:ny-1, 0:nx-1, :))

            call write_h5_dataset(h5file_id, name, "total buoyancy", &
                                  tbuoyg(0:nz, 0:ny-1, 0:nx-1))

#ifndef ENABLE_DRY_MODE
            call write_h5_dataset(h5file_id, name, "dry buoyancy", &
                                  dbuoyg(0:nz, 0:ny-1, 0:nx-1))
#endif
#ifdef ENABLE_DIAGNOSE
            call write_h5_dataset(h5file_id, name, "volume", &
                                  volg(0:nz, 0:ny-1, 0:nx-1))

            call write_h5_dataset(h5file_id, name, "nparg", &
                                  nparg(0:nz-1, :, :))
#endif

#ifndef NDEBUG
            call write_h5_dataset(h5file_id, name, "symmetry volume", &
                                  sym_volg(0:nz, 0:ny-1, 0:nx-1))

            call write_h5_dataset(h5file_id, name, "velocity gradient tensor", &
                                  velgradg(0:nz, 0:ny-1, 0:nx-1, :))

            call write_h5_dataset(h5file_id, name, "vorticity tendency", &
                                  vtend(0:nz, 0:ny-1, 0:nx-1, :))
#endif

            call close_h5_group(group)
        end subroutine write_h5_fields

        subroutine read_h5_domain(h5fname)
            character(*), intent(in) :: h5fname

            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5file_id)
            call read_h5_box(h5file_id, ncells, extent, lower)
            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)
            call close_h5_file(h5file_id)

        end subroutine read_h5_domain

end module field_hdf5
