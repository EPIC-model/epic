module fields
    use parameters, only : mesh
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_3d, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none

    double precision, allocatable, dimension(:, :, :) :: &
        velocity_f,       &   ! velocity vector field
        strain_f,         &   ! velocity gradient tensor
        volume_f              ! volume scalar field

    contains
        function get_mesh_spacing() result(dx)
            double precision :: dx(2)
            dx = mesh%extent / (mesh%grid - 1)
        end function get_mesh_spacing

        ! get the lower field index given the parcel position
        function get_lower_index(pos) result(idx)
            double precision, intent(in)  :: pos(2)
            integer                       :: idx(2)
            double precision              :: dx(2)

            dx = get_mesh_spacing()

            ! + 1 since Fortran starts with 1
            idx = (pos - mesh%origin) / dx + 1

        end function get_lower_index

        ! get a position given a field index
        function get_position(idx) result(pos)
            integer,         intent(in) :: idx(2)
            double precision            :: pos(2)
            double precision            :: dx(2)

            dx = get_mesh_spacing()

            ! we need to subtract 1 from the index
            ! since Fortran starts with 1 not with 0
            pos = mesh%origin + (idx - 1) * dx

        end function get_position


        subroutine write_h5_fields(iter)
            integer, intent(in)        :: iter ! iteration
            integer(hid_t)             :: group
            integer(hid_t)             :: step_group
            character(:), allocatable  :: step
            character(:), allocatable  :: name
            logical                    :: link_exists = .false.

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/fields"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            !
            ! write fields
            !

            call write_h5_dataset_3d(name, "velocity", velocity_f)

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_fields

end module fields
