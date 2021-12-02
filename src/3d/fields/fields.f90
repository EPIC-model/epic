! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use parameters, only : dx, dxi, extent, lower, nx, ny, nz
    use constants, only : zero
    implicit none

    ! x: zonal
    ! y: meridional
    ! z: vertical
    ! Halo grid points in vertical direction z are -1 and nz+1,
    ! hence the valid regrion is from 0 to nz
    ! Due to periodicity in x and y, the grid points in x go from 0 to nx-1
    ! and from 0 to ny-1 in y
    double precision, allocatable, dimension(:, :, :, :) :: &
        velog,     &   ! velocity vector field
        vortg,     &   ! vorticity vector field
        vtend,     &   ! vorticity tendency
        velgradg       ! velocity gradient tensor
                       ! ordering: du/dx, du/dy,
                       !                  dv/dy,
                       !           dw/dx, dw/dy
                       ! the derivatives dv/dx, du/dz, dv/dz and dw/dz
                       ! are calculated on the fly with vorticity
                       ! or the assumption of incompressibility (du/dx + dv/dy + dw/dz = 0):
                       !    dv/dx = \omegaz + du/dy
                       !    du/dz = \omegay + dw/dx
                       !    dv/dz = dw/dy - \omegax
                       !    dw/dz = - (du/dx + dv/dy)

    double precision, allocatable, dimension(:, :, :) :: &
#ifndef ENABLE_DRY_MODE
        dbuoyg,    &   ! dry buoyancy (or liquid-water buoyancy)
#endif
        tbuoyg,    &   ! buoyancy
#ifndef NDEBUG
        sym_volg,  &   ! symmetry volume (debug mode only)
#endif
        volg           ! volume scalar field

    integer, allocatable, dimension(:, :, :) :: &
        nparg,     &   ! number of parcels per grid box
        nsparg         ! number of small parcels per grid box

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:ny-1, 0:nx-1, 3))
            allocate(velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5))

            allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

#ifndef NDEBUG
            allocate(sym_volg(-1:nz+1, 0:ny-1, 0:nx-1))
#endif

            allocate(vortg(-1:nz+1, 0:ny-1, 0:nx-1, 3))

            allocate(vtend(-1:nz+1, 0:ny-1, 0:nx-1, 3))

            allocate(tbuoyg(-1:nz+1, 0:ny-1, 0:nx-1))

#ifndef ENABLE_DRY_MODE
            allocate(dbuoyg(-1:nz+1, 0:ny-1, 0:nx-1))
#endif

            allocate(nparg(-1:nz, 0:ny-1, 0:nx-1))
            allocate(nsparg(-1:nz, 0:ny-1, 0:nx-1))

        end subroutine field_alloc

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            volg     = zero
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
        ! @param[out] i lower, zonal cell index
        ! @param[out] j lower, vertical cell index
        pure subroutine get_index(pos, i, j, k)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: i, j, k

            i = floor((pos(1) - lower(1)) * dxi(1))
            j = floor((pos(2) - lower(2)) * dxi(2))
            k = floor((pos(3) - lower(3)) * dxi(3))
        end subroutine get_index


        ! Do periodic shift of the index
        ! @param[inout] ii zonal grid point indices
        ! @param[inout] jj meridional grid point indices
        pure subroutine periodic_index_shift(ii, jj)
            integer, intent(inout) :: ii(:), jj(:)

            ! account for x / y periodicity:
            ! -1          --> nx-1 / ny-1
            !  0          --> 0
            ! nx+1 / ny+1 --> 1
            ! nx / ny     --> 0
            ! nx-1 / ny-1 --> nx-1 / ny-1
            ii = mod(ii + nx, nx)
            jj = mod(jj + ny, ny)

        end subroutine periodic_index_shift


        ! Get the coordinate of a grid point (i, j, k).
        ! @param[in] i zonal cell index
        ! @param[in] j meridional cell index
        ! @param[in] k vertical cell index
        ! @param[out] pos position of (i, j, k) in the domain
        pure subroutine get_position(i, j, k, pos)
            integer,          intent(in)  :: i, j, k
            double precision, intent(out) :: pos(3)

            pos(1) = lower(1) + i * dx(1)
            pos(2) = lower(2) + j * dx(2)
            pos(3) = lower(3) + k * dx(3)

        end subroutine get_position

end module fields
