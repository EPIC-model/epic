! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use dimensions, only : n_dim, I_X, I_Y, I_Z
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
        velog,     &   ! velocity vector field (u, v, w)
        vortg,     &   ! vorticity vector field (\xi, \eta, \zeta)
        vtend,     &   ! vorticity tendency
        velgradg       ! velocity gradient tensor
                       ! ordering: du/dx, du/dy,
                       !                  dv/dy,
                       !           dw/dx, dw/dy
                       ! the derivatives dv/dx, du/dz, dv/dz and dw/dz
                       ! are calculated on the fly with vorticity
                       ! or the assumption of incompressibility (du/dx + dv/dy + dw/dz = 0):
                       !    dv/dx = \zeta + du/dy
                       !    du/dz = \eta + dw/dx
                       !    dv/dz = dw/dy - \xi
                       !    dw/dz = - (du/dx + dv/dy)

    double precision, allocatable, dimension(:, :, :) :: &
#ifndef ENABLE_DRY_MODE
        dbuoyg,    &   ! dry buoyancy (or liquid-water buoyancy)
        humg,      &   ! humidity
#endif
        tbuoyg,    &   ! buoyancy
#ifndef NDEBUG
        sym_volg,  &   ! symmetry volume (debug mode only)
#endif
        volg           ! volume scalar field

    integer, allocatable, dimension(:, :, :) :: &
        nparg,     &   ! number of parcels per grid box
        nsparg         ! number of small parcels per grid box

    ! velocity strain indices
    integer, parameter :: I_DUDX = 1 & ! index for du/dx strain component
                        , I_DUDY = 2 & ! index for du/dy strain component
                        , I_DVDY = 3 & ! index for dv/dy strain component
                        , I_DWDX = 4 & ! index for dw/dx strain component
                        , I_DWDY = 5   ! index for dw/dy strain component

    contains

        ! Allocate all fields
        subroutine field_alloc
            if (allocated(velog)) then
                return
            endif

            allocate(velog(-1:nz+1, 0:ny-1, 0:nx-1, n_dim))
            allocate(velgradg(-1:nz+1, 0:ny-1, 0:nx-1, 5))

            allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

#ifndef NDEBUG
            allocate(sym_volg(-1:nz+1, 0:ny-1, 0:nx-1))
#endif

            allocate(vortg(-1:nz+1, 0:ny-1, 0:nx-1, n_dim))

            allocate(vtend(-1:nz+1, 0:ny-1, 0:nx-1, n_dim))

            allocate(tbuoyg(-1:nz+1, 0:ny-1, 0:nx-1))

#ifndef ENABLE_DRY_MODE
            allocate(dbuoyg(-1:nz+1, 0:ny-1, 0:nx-1))
            allocate(humg(-1:nz+1, 0:ny-1, 0:nx-1))
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
            humg     = zero
#endif
            nparg    = zero
            nsparg   = zero
        end subroutine

        ! Get the lower index of the cell the parcel is in.
        ! This subroutine does not take x periodicity into account.
        ! @param[in] pos position of the parcel
        ! @param[out] i lower, zonal cell index
        ! @param[out] j lower, meridional cell index
        ! @param[out] k lower, vertical cell index
        pure subroutine get_index(pos, i, j, k)
            double precision, intent(in)  :: pos(n_dim)
            integer,          intent(out) :: i, j, k

            i = floor((pos(I_X) - lower(I_X)) * dxi(I_X))
            j = floor((pos(I_Y) - lower(I_Y)) * dxi(I_Y))
            k = floor((pos(I_Z) - lower(I_Z)) * dxi(I_Z))
        end subroutine get_index

        ! Get the lower horizontal (x, y) index of the cell the parcel is in.
        ! This subroutine does not take x periodicity into account.
        ! @param[in] pos position of the parcel
        ! @param[out] i lower, zonal cell index
        ! @param[out] j lower, meridional cell index
        pure subroutine get_horizontal_index(pos, i, j)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: i, j

            i = floor((pos(I_X) - lower(I_X)) * dxi(I_X))
            j = floor((pos(I_Y) - lower(I_Y)) * dxi(I_Y))
        end subroutine get_horizontal_index


        ! Do periodic shift of the index
        ! @param[inout] ii zonal grid point indices
        ! @param[inout] jj meridional grid point indices
        elemental pure subroutine periodic_index_shift(ii, jj)
            integer, intent(inout) :: ii, jj

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
            double precision, intent(out) :: pos(n_dim)

            pos(I_X) = lower(I_X) + i * dx(I_X)
            pos(I_Y) = lower(I_Y) + j * dx(I_Y)
            pos(I_Z) = lower(I_Z) + k * dx(I_Z)

        end subroutine get_position

        ! Get the coordinate of a grid point (i, j).
        ! @param[in] i zonal cell index
        ! @param[in] j meridional cell index
        ! @param[out] pos position of (i, j) in the domain
        pure subroutine get_horizontal_position(i, j, pos)
            integer,          intent(in)  :: i, j
            double precision, intent(out) :: pos(2)

            pos(I_X) = lower(I_X) + i * dx(I_X)
            pos(I_Y) = lower(I_Y) + j * dx(I_Y)

        end subroutine get_horizontal_position

end module fields
