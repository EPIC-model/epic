! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_ellipse_interpl
    use constants, only : zero, one, two, f12
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, ny, amin, ncelli
    use options, only : parcel
    use fields
    use parcels_mod
    use omp_lib
    implicit none

    private

    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer :: is(ngp), js(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

!     integer :: surf_par2grid_timer, &
!                surf_grid2par_timer

!     private :: is, js, weights

    public :: area2grid

    contains

        subroutine area2grid
            call m_area2grid(0,  bot_parcels)
            call m_area2grid(nz, top_parcels)
        end subroutine area2grid

        ! Interpolate the parcel area to the grid
        subroutine m_area2grid(iz, surf_parcels)
            integer,               intent(in) :: iz
            type(ellipse_pc_type), intent(in) :: surf_parcels
            double precision                  :: points(2, 2)
            integer                           :: n, p, l

            volg(iz, :, :) = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, is, js, weights) &
            !$omp& reduction(+: volg)
            do n = 1, surf_parcels%local_num
                points = surf_parcels%get_points(n)

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
!                     call apply_surface_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
!                     call bilinear(points(:, p), is, js, weights)

                    do l = 1, ngp
                        volg(iz, js(l), is(l)) = volg(iz, js(l), is(l)) &
                                              + f12 * weights(l) * surf_parcels%area(n)
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine m_area2grid

!         subroutine surface_par2grid
!             call do_par2grid(lo_surf_parcels, n_lo_surf_parcels, 'lo')
!             call do_par2grid(up_surf_parcels, n_up_surf_parcels, 'up')
!         end subroutine surface_par2grid
!
!         ! Interpolate parcel quantities to the grid, these consist of the parcel
!         !   - potential vorticity
!         !   - area
!         !   - divergence
!         !   - dimensionless height
!         ! It also updates the scalar fields:
!         !   - nparg, that is the number of parcels per grid cell
!         !   - nsparg, that is the number of small parcels per grid cell
!         subroutine do_par2grid(s_parcels, n_par, which)
!             type(surface_parcel_container_type), intent(in) :: s_parcels
!             integer,                             intent(in) :: n_par
!             character(2),                        intent(in) :: which
!             double precision                                :: points(2, 2)
!             integer                                         :: n, p, l, k
!             double precision                                :: weight, btot
!
!             call start_timer(surf_par2grid_timer)
!
!             k = nz
!             if (which == 'lo') then
!                 k = 0
!             endif
!
!             volg(k, :, :) = zero
!             vortg(k, :, :, :) = zero
!             tbuoyg(k, :, :) = zero
!
! !             nparg = zero
! !             nsparg = zero
!
!             !$omp parallel default(shared)
!             !$omp do private(n, p, l, points, weight, is, js, weights, btot) &
!             !$omp& reduction(+:volg, vortg, tbuoyg)
!             do n = 1, n_par
!                 points = get_ellipse_points(s_parcels%position(:, n), &
!                                             s_parcels%B(:, n))
!
! !                 call get_index(s_parcels%position(:, n), i, j)
! !                 i = mod(i + nx, nx)
! !                 j = mod(j + ny, ny)
! !                 nparg(j, i) = nparg(j, i) + 1
! !                 if (s_parcels%area(n) <= amin) then
! !                     nsparg(j, i) = nsparg(j, i) + 1
! !                 endif
!
!                 btot = s_parcels%buoyancy(n)
!
!                 ! we have 2 points per ellipse
!                 do p = 1, 2
!
!                     ! ensure point is within the domain
!                     call apply_surface_periodic_bc(points(:, p))
!
!                     ! get interpolation weights and mesh indices
!                     call bilinear(points(:, p), is, js, weights)
!
!                     ! loop over grid points which are part of the interpolation
!                     ! the weight is halved due to 2 points per ellipse
!                     do l = 1, ngp
!
!                         weight = f12 * weights(l) * s_parcels%area(n)
!
!                         volg(k, js(l), is(l)) = volg(k, js(l), is(l)) + weight
!
!                         vortg(k, js(l), is(l), :) = vortg(k, js(l), is(l), :) &
!                                                   + weight * s_parcels%vorticity(:, n)
!
!                         tbuoyg(k, js(l), is(l)) = tbuoyg(k, js(l), is(l)) &
!                                                 + weight * btot
!
!                     enddo
!                 enddo
!             enddo
!             !$omp end do
!             !$omp end parallel
!
! !             ! sanity check
! !             if (sum(nparg) /= n_par) then
! !                 print *, "par2grid: Wrong total number of parcels!"
! !                 stop
! !             endif
!
!             call stop_timer(surf_par2grid_timer)
!
!         end subroutine do_par2grid
!
!
!         ! Interpolate the gridded quantities to the parcels
!         ! @param[inout] vel is the parcel velocity
!         ! @param[inout] vgrad is the parcel strain
!         ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
!         !            (optional)
!         subroutine do_grid2par(s_parcels, n_par, which, vel, vor, vgrad, add)
!             type(surface_parcel_container_type), intent(inout) :: s_parcels
!             integer,                             intent(in)    :: n_par
!             character(2),                        intent(in)    :: which
!             double precision,     intent(inout) :: vel(:, :), vor(:, :), vgrad(:, :)
!             logical, optional, intent(in)       :: add
!             double precision                    :: points(2, 2), weight, dvdx
!             integer                             :: n, p, l, k
!
!             call start_timer(surf_grid2par_timer)
!
!             k = nz
!             if (which == 'lo') then
!                 k = 0
!             endif
!
!             ! clear old data efficiently
!             if(present(add)) then
!                if(add .eqv. .false.) then
!                     !$omp parallel default(shared)
!                     !$omp do private(n)
!                     do n = 1, n_par
!                         vel(:, n) = zero
!                         vor(:, n) = zero
!                     enddo
!                     !$omp end do
!                     !$omp end parallel
!                endif
!             else
!                 !$omp parallel default(shared)
!                 !$omp do private(n)
!                 do n = 1, n_par
!                     vel(:, n) = zero
!                     vor(:, n) = zero
!                 enddo
!                 !$omp end do
!                 !$omp end parallel
!             endif
!
!             !$omp parallel default(shared)
!             !$omp do private(n, p, l, points, weight, is, js, weights, dvdx)
!             do n = 1, n_par
!
!                 vgrad(:, n) = zero
!
!                 points = get_ellipse_points(s_parcels%position(:, n), &
!                                             s_parcels%B(:, n))
!
!                 ! we have 2 points per ellipse
!                 do p = 1, 2
!
!                     ! ensure point is within the domain
!                     call apply_surface_periodic_bc(points(:, p))
!
!                     ! get interpolation weights and mesh indices
!                     call bilinear(points(:, p), is, js, weights)
!
!                     ! loop over grid points which are part of the interpolation
!                     do l = 1, ngp
!                         ! the weight is halved due to 2 points per ellipse
!                         weight = f12 * weights(l)
!
!                         vel(:, n) = vel(:, n) + weight * velog(k, js(l), is(l), 1:2)
!
!                         ! du/dx
!                         vgrad(1, n) = vgrad(1, n) + weight * velgradg(k, js(l), is(l), 1)
!
!                         ! du/dy
!                         vgrad(2, n) = vgrad(2, n) + weight * velgradg(k, js(l), is(l), 2)
!
!                         ! dv/dx = \zeta + du/dy
!                         dvdx = vortg(k, js(l), is(l), 3) + velgradg(k, js(l), is(l), 2)
!                         vgrad(3, n) = vgrad(3, n) + weight * dvdx
!
!                         ! dv/dy
!                         vgrad(4, n) = vgrad(4, n) + weight * velgradg(k, js(l), is(l), 3)
!
!                         vor(:, n) = vor(:, n) + weight * vtend(k, js(l), is(l), :)
!
!                     enddo
!                 enddo
!             enddo
!             !$omp end do
!             !$omp end parallel
!
!             call stop_timer(surf_grid2par_timer)
!
!         end subroutine do_grid2par

!         ! Bi-linear interpolation
!         ! @param[in] pos position of the parcel
!         ! @param[out] ii horizontal grid points for interoplation
!         ! @param[out] jj vertical grid points for interpolation
!         ! @param[out] ww interpolation weights
!         subroutine bilinear(pos, ii, jj, ww)
!             double precision, intent(in)  :: pos(2)
!             integer,          intent(out) :: ii(4), jj(4)
!             double precision, intent(out) :: ww(4)
!             double precision              :: xy(2)
!
!             ! (i, j)
!             call get_surf_index(pos, ii(1), jj(1))
!             call get_surf_position(ii(1), jj(1), xy)
!             ww(1) = product(one - abs(pos - xy) * dxi(1:2))
!
!             ! (i+1, j)
!             ii(2) = ii(1) + 1
!             jj(2) = jj(1)
!             call get_surf_position(ii(2), jj(2), xy)
!             ww(2) = product(one - abs(pos - xy) * dxi(1:2))
!
!             ! (i, j+1)
!             ii(3) = ii(1)
!             jj(3) = jj(1) + 1
!             call get_surf_position(ii(3), jj(3), xy)
!             ww(3) = product(one - abs(pos - xy) * dxi(1:2))
!
!             ! (i+1, j+1)
!             ii(4) = ii(2)
!             jj(4) = jj(3)
!             call get_surf_position(ii(4), jj(4), xy)
!             ww(4) = product(one - abs(pos - xy) * dxi(1:2))
!
!             ! account for x periodicity
!             call periodic_index_shift(ii, jj)
!
!         end subroutine bilinear

end module parcel_ellipse_interpl
