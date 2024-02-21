! =============================================================================
!     This module specifies all fields and implements specific subroutines
!     and functions.
! =============================================================================
module fields
    use dimensions, only : n_dim, I_X, I_Y, I_Z
    use parameters, only : dx, dxi, extent, lower, nx, ny, nz
    use constants, only : zero
    use mpi_environment
    use mpi_layout, only : box, l_mpi_layout_initialised
    use mpi_utils, only : mpi_exit_on_error
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
                       ! ordering: du/dx, du/dy, du/dz
                       !           dv/dx, dv/dy, dv/dz
                       !           dw/dx, dw/dy
                       ! the derivative dw/dz
                       ! is calculated with
                       ! the assumption of incompressibility (du/dx + dv/dy + dw/dz = 0):
                       !    dw/dz = - (du/dx + dv/dy)
!                        !    dv/dx = \zeta + du/dy
!                        !    du/dz = \eta + dw/dx
!                        !    dv/dz = dw/dy - \xi

    double precision, allocatable, dimension(:, :, :) :: &
#ifndef ENABLE_DRY_MODE
        dbuoyg,    &   ! dry buoyancy (or liquid-water buoyancy)
        humg,      &   ! humidity
#endif
        tbuoyg,    &   ! buoyancy
#ifndef NDEBUG
        sym_volg,  &   ! symmetry volume (debug mode only)
#endif
        volg, &        ! volume scalar field
        strain_mag     ! strain magnitude

    integer, allocatable, dimension(:, :, :) :: &
        nparg,     &   ! number of parcels per grid box
        nsparg         ! number of small parcels per grid box

    ! velocity strain indices (note that dw/dz is found from dw/dz = - (du/dx + dv/dy))
    integer, parameter :: I_DUDX = 1 & ! index for du/dx strain component
                        , I_DUDY = 2 & ! index for du/dy strain component
                        , I_DUDZ = 3 & ! index for du/dz strain component
                        , I_DVDX = 4 & ! index for dv/dx strain component
                        , I_DVDY = 5 & ! index for dv/dy strain component
                        , I_DVDZ = 6 & ! index for dv/dz strain component
                        , I_DWDX = 7 & ! index for dw/dx strain component
                        , I_DWDY = 8   ! index for dw/dy strain component


    contains

        ! Allocate all fields
        subroutine field_alloc
            integer :: hlo(3), hhi(3)

            if (.not. l_mpi_layout_initialised) then
                call mpi_exit_on_error
            endif

            if (allocated(velog)) then
                return
            endif

            hlo = box%hlo
            hhi = box%hhi

            allocate(velog(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1), n_dim))
            allocate(velgradg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1), 8))

            allocate(volg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))
            allocate(strain_mag(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))

#ifndef NDEBUG
            allocate(sym_volg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))
#endif

            allocate(vortg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1), n_dim))

            allocate(vtend(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1), n_dim))

            allocate(tbuoyg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))

#ifndef ENABLE_DRY_MODE
            allocate(dbuoyg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))
            allocate(humg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))
#endif

            allocate(nparg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))
            allocate(nsparg(hlo(3):hhi(3), hlo(2):hhi(2), hlo(1):hhi(1)))

        end subroutine field_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Reset fields to zero
        subroutine field_default
            call field_alloc

            velog    = zero
            velgradg = zero
            volg     = zero
            strain_mag = zero
            vortg    = zero
            vtend    = zero
            tbuoyg   = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg   = zero
            humg     = zero
#endif
            nparg    = zero
            nsparg   = zero

#ifndef NDEBUG
            sym_volg = zero
#endif
        end subroutine field_default

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate fields
        subroutine field_dealloc
            if (allocated(velog)) then
                deallocate(velog)
                deallocate(velgradg)
                deallocate(volg)
                deallocate(strain_mag)
                deallocate(vortg)
                deallocate(vtend)
                deallocate(tbuoyg)
#ifndef ENABLE_DRY_MODE
                deallocate(dbuoyg)
                deallocate(humg)
#endif
                deallocate(nparg)
                deallocate(nsparg)

#ifndef NDEBUG
                deallocate(sym_volg)
#endif
            endif

        end subroutine field_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Get the lower index of the cell the parcel is in.
        ! This subroutine does not take x nor y periodicity into account.
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        pure function is_contained(pos) result(l_contained)
            double precision, intent(in) :: pos(3)
            integer                      :: i, j, k
            logical                      :: l_contained

            call get_index(pos, i, j, k)

            l_contained = ((i >= box%lo(1))  .and. &
                           (i <= box%hi(1))  .and. &
                           (j >= box%lo(2))  .and. &
                           (j <= box%hi(2)))
        end function

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

end module fields
