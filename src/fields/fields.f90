! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : dx, dxi, extent, lower, nx, nz
    use constants, only : zero
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_2d, &
                       write_h5_dataset_3d, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none

    ! Halo grid points in vertical direction z are -1 and nz+1,
    ! hence the valid regrion is from 0 to nz
    ! Due to periodicity in x, the grid points in x go from 0 to nx-1
    double precision, allocatable, dimension(:, :, :) :: &
        velog,     &   ! velocity vector field (has 1 halo cell layer in z)
        velgradg,  &   ! velocity gradient tensor (has 1 halo cell layer in z)
        buoyg,     &   ! buoyancy (has 1 halo cell layer in z)
        humg,      &   ! specific humidity
        humlig,    &   ! condensed humidity
        vortg          ! vorticity scalar field (has no halo cell layers)

    double precision, allocatable, dimension(:, :) :: &
        volg           ! volume scalar field (has 1 halo cell layer in z)

    contains

        ! allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:nx-1, 2))
            allocate(velgradg(-1:nz+1, 0:nx-1, 4))

            allocate(volg(-1:nz+1, 0:nx-1))

            ! vorticity has no halo grid points in y
            allocate(vortg(-1:nz+1, 0:nx-1, 1))

            allocate(buoyg(-1:nz+1, 0:nx-1, 1))

            allocate(humg(-1:nz+1, 0:nx-1, 1))

            allocate(humlig(-1:nz+1, 0:nx-1, 1))

        end subroutine field_alloc

        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            volg     = zero
            buoyg    = zero
            humg     = zero
            humlig   = zero
        end subroutine

        ! get the lower index of the cell the parcel is in
        ! this subroutine does not take x periodicity into account
        subroutine get_index(pos, i, j)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: i, j
            integer                       :: idx(2)

            idx = floor((pos - lower) * dxi)

            i = idx(1)
            j = idx(2)
        end subroutine get_index


        ! do periodic shift of the index
        subroutine periodic_index_shift(ii)
            integer, intent(inout) :: ii(:)

            ! account for x periodicity:
            ! -1   --> nx-1
            !  0   --> 0
            ! nx+1 --> 1
            ! nx   --> 0
            ! nx-1 --> nx-1
            ii = mod(ii + nx, nx)

        end subroutine periodic_index_shift


        ! get a position given a field index
        subroutine get_position(i, j, pos)
            integer,          intent(in)  :: i, j
            double precision, intent(out) :: pos(2)

            pos = lower + (/dble(i), dble(j)/) * dx

        end subroutine get_position


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
            ! write fields (do not write halo cells)
            !

            call write_h5_dataset_3d(name, "velocity", &
                                     velog(0:nz, 0:nx-1, :))

            call write_h5_dataset_3d(name, "velocity gradient tensor", &
                                     velgradg(0:nz, 0:nx-1, :))

            call write_h5_dataset_2d(name, "volume", &
                                     volg(0:nz, 0:nx-1))

            call write_h5_dataset_3d(name, "buoyancy", &
                                     buoyg(0:nz, 0:nx-1, :))

            call write_h5_dataset_3d(name, "humidity", &
                                     humg(0:nz, 0:nx-1, :))

            call write_h5_dataset_3d(name, "vorticity", &
                                     vortg(0:nz, 0:nx-1, :))

            call write_h5_dataset_3d(name, "liquid humidity", &
                                     humlig(0:nz, 0:nx-1, :))

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_fields

end module fields
