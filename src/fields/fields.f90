module fields
    use options, only : grid
    use parameters, only : dx, dxi, extent, lower
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_3d, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none

    ! Halo grid points are -1 and grid(1) and grid(2),
    ! hence the valid regrion is from 0 to grid(1) or grid(2), respectively.
    double precision, allocatable, dimension(:, :, :) :: &
        velocity_f,       &   ! velocity vector field (has 1 halo cell layer)
        strain_f,         &   ! velocity gradient tensor (has 1 halo cell layer)
        volume_f,         &   ! volume scalar field (has 1 halo cell layer)
        vorticity_f           ! vorticity scalar field (has no halo cell layers)

    contains

        ! get the lower index of the cell the parcel is in
        ! this function does not take x periodicity into account
        function get_index(pos) result(idx)
            double precision, intent(in)  :: pos(2)
            integer                       :: idx(2)

            idx = (pos - lower) * dxi

        end function get_index


        ! do periodic shift of the index
        subroutine periodic_index_shift(idx, n)
            integer, intent(inout) :: idx(2, n)
            integer, intent(in)    :: n

            ! account for x periodicity:
            ! [nx = grid(1)]
            ! -1   --> nx-2
            !  0   --> 0
            ! nx   --> 1
            ! nx-1 --> 0
            idx(1, :) = mod(idx(1, :) + grid(1) - 1, grid(1) - 1)

        end subroutine periodic_index_shift


        ! get a position given a field index
        function get_position(idx) result(pos)
            integer,         intent(in) :: idx(2)
            double precision            :: pos(2)

            pos = lower + idx * dx

        end function get_position


        subroutine write_h5_fields(iter)
            integer, intent(in)        :: iter ! iteration
            integer(hid_t)             :: group
            integer(hid_t)             :: step_group
            character(:), allocatable  :: step
            character(:), allocatable  :: name

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/fields"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            !
            ! write fields
            !

            if (iter == 0) then
                ! do not write halo cells
                call write_h5_dataset_3d(name, "velocity",          &
                    velocity_f(0:grid(1)-1, 0:grid(2)-1, :))

                ! do not write halo cells
                call write_h5_dataset_3d(name, "velocity strain",   &
                    strain_f(0:grid(1)-1, 0:grid(2)-1, :))
            endif

            ! do not write halo cells
            call write_h5_dataset_3d(name, "volume",            &
                volume_f(0:grid(1)-1, 0:grid(2)-1, :))

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_fields

end module fields
