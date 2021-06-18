! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : dx, dxi, extent, lower, nx, nz
    use constants, only : zero
    implicit none

    ! Halo grid points in vertical direction z are -1 and nz+1,
    ! hence the valid regrion is from 0 to nz
    ! Due to periodicity in x, the grid points in x go from 0 to nx-1
    double precision, allocatable, dimension(:, :, :) :: &
        velog,     &   ! velocity vector field (has 1 halo cell layer in z)
        velgradg       ! velocity gradient tensor (has 1 halo cell layer in z)

    double precision, allocatable, dimension(:, :) :: &
        vortg,     &   ! vorticity scalar field
        vtend,     &   ! vorticity tendency
        tbuoyg,    &   ! buoyancy (has 1 halo cell layer in z)
        humg,      &   ! specific humidity
        humlig,    &   ! condensed humidity
#ifndef NDEBUG
        sym_volg,  &   ! symmetry volume (debug mode only)
#endif
        volg           ! volume scalar field (has 1 halo cell layer in z)

    integer, allocatable, dimension(:, :) :: &
        nparg          ! number of parcels per grid box (from 0 to nz-1 and 0 to nx-1)

    contains

        ! allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:nx-1, 2))
            allocate(velgradg(-1:nz+1, 0:nx-1, 4))

            allocate(volg(-1:nz+1, 0:nx-1))

#ifndef NDEBUG
            allocate(sym_volg(-1:nz+1, 0:nx-1))
#endif

            allocate(vortg(-1:nz+1, 0:nx-1))

            allocate(vtend(-1:nz+1, 0:nx-1))

            allocate(tbuoyg(-1:nz+1, 0:nx-1))

            allocate(humg(-1:nz+1, 0:nx-1))

            allocate(humlig(-1:nz+1, 0:nx-1))

            allocate(nparg(-1:nz, 0:nx-1))

        end subroutine field_alloc

        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            volg     = zero
            vortg    = zero
            vtend    = zero
            tbuoyg    = zero
            humg     = zero
            humlig   = zero
            nparg    = zero
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

end module fields
