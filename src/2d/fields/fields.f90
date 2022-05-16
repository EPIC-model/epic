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
        velog,     &   ! velocity vector field
        velgradg       ! velocity gradient tensor

    double precision, allocatable, dimension(:, :) :: &
        vortg,     &   ! vorticity scalar field
        vtend,     &   ! vorticity tendency
#ifndef ENABLE_DRY_MODE
        dbuoyg,    &   ! dry buoyancy (or liquid-water buoyancy)
#endif
        tbuoyg,    &   ! buoyancy
#ifndef NDEBUG
        sym_volg,  &   ! symmetry volume (debug mode only)
#endif
        volg,      &   ! volume scalar field
        humg

    integer, allocatable, dimension(:, :) :: &
        nparg,     &   ! number of parcels per grid box
        nsparg         ! number of small parcels per grid box

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:nx-1, 2))
            allocate(velgradg(-1:nz+1, 0:nx-1, 4))

            allocate(volg(-1:nz+1, 0:nx-1))

            allocate(humg(-1:nz+1, 0:nx-1))

#ifndef NDEBUG
            allocate(sym_volg(-1:nz+1, 0:nx-1))
#endif

            allocate(vortg(-1:nz+1, 0:nx-1))

            allocate(vtend(-1:nz+1, 0:nx-1))

            allocate(tbuoyg(-1:nz+1, 0:nx-1))

#ifndef ENABLE_DRY_MODE
            allocate(dbuoyg(-1:nz+1, 0:nx-1))
#endif

            allocate(nparg(-1:nz, 0:nx-1))
            allocate(nsparg(-1:nz, 0:nx-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            volg     = zero
            humg      = zero
            vortg    = zero
            vtend    = zero
            tbuoyg   = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg   = zero
#endif
            nparg    = zero
            nsparg   = zero
        end subroutine

        ! Get the lower index of the cell the parcel is in.
        ! This subroutine does not take x periodicity into account.
        ! @param[in] pos position of the parcel
        ! @param[out] i lower, horizontal cell index
        ! @param[out] j lower, vertical cell index
        subroutine get_index(pos, i, j)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: i, j
            integer                       :: idx(2)

            idx = floor((pos - lower) * dxi)

            i = idx(1)
            j = idx(2)
        end subroutine get_index


        ! Do periodic shift of the index
        ! @param[inout] ii horizontal grid point indices
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


        ! Get the coordinate of a grid point (i, j).
        ! @param[in] i horizontal cell index
        ! @param[in] j vertical cell index
        ! @param[out] pos position of (i, j) in the domain
        subroutine get_position(i, j, pos)
            integer,          intent(in)  :: i, j
            double precision, intent(out) :: pos(2)

            pos = lower + (/dble(i), dble(j)/) * dx

        end subroutine get_position

end module fields
