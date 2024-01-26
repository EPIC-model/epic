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
    use interpl, only : bilinear
    use parcel_bc, only : apply_periodic_bc
    use physics, only : glat, lambda_c, q_0
#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
    use physics, only : bfsq
#endif

    use omp_lib
    implicit none

    private

!     integer :: surf_par2grid_timer, &
!                surf_grid2par_timer

!     private :: is, js, weights

#ifndef ENABLE_G2P_1POINT
    integer, parameter :: n_points_g2p = 2
    double precision, parameter :: point_weight_g2p = f12
#else
    integer, parameter :: n_points_g2p = 1
    double precision, parameter :: point_weight_g2p = one
#endif

#ifndef ENABLE_P2G_1POINT
    integer, parameter :: n_points_p2g = 2
    double precision, parameter :: point_weight_p2g = f12
#else
    integer, parameter :: n_points_p2g = 1
    double precision, parameter :: point_weight_p2g = one
#endif


    public :: area2grid         &
            , surf_par2grid     &
            , surf_grid2par

    contains

        subroutine area2grid
            call m_area2grid(0,  bot_parcels)
            call m_area2grid(nz, top_parcels)
        end subroutine area2grid

        ! Interpolate the parcel area to the grid
        ! ATTENTION After this operation the halo grid points
        !           are not filled properly.
        subroutine m_area2grid(iz, surf_parcels)
            integer,               intent(in) :: iz
            type(ellipse_pc_type), intent(in) :: surf_parcels
            double precision                  :: points(2, 2)
            integer                           :: n, p, l, is, js
            double precision                  :: weights(0:1, 0:1)

            volg(iz, :, :) = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, is, js, weights) &
            !$omp& reduction(+: volg)
            do n = 1, surf_parcels%local_num
                points = surf_parcels%get_points(n)

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    volg(iz, js:js+1, is:is+1) = volg(iz, js:js+1, is:is+1) &
                                               + f12 * weights * surf_parcels%area(n)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine m_area2grid

        ! ATTENTION After this operation the halo grid points
        !           are not filled properly.
        subroutine surf_par2grid
!             call start_timer(surf_par2grid_timer)
            call m_par2grid(0,  bot_parcels)
            call m_par2grid(nz, top_parcels)
!             call stop_timer(surf_par2grid_timer)
        end subroutine surf_par2grid

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !  - vorticity
        !  - area
        !  - divergence
        ! It also updates the scalar fields:
        !  - nparg, that is the number of parcels per grid cell
        !  - nsparg, that is the number of small parcels per grid cell
        ! ATTENTION After this operation the halo grid points
        !           are not filled properly.
        subroutine m_par2grid(iz, surf_parcels)
            integer,               intent(in) :: iz
            type(ellipse_pc_type), intent(in) :: surf_parcels
            double precision                  :: points(2, n_points_p2g)
            integer                           :: n, p, l, is, js
            double precision                  :: weights(0:1, 0:1), btot
            double precision                  :: zpos
#ifndef ENABLE_DRY_MODE
            double precision                  :: qq, q_c
#endif

            vortg(iz, :, :, :) = zero
            volg(iz, :, :) = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg(iz, :, :) = zero
            humg(iz, :, :) = zero
#endif
            tbuoyg(iz, :, :) = zero

            nparg = zero
            nsparg = zero

            zpos = lower(3) + dble(iz) * dx(3)
#ifndef ENABLE_DRY_MODE
            qq = q_0 * dexp(lambda_c * (lower(3) - zpos))
#endif

            !$omp parallel default(shared)
#ifndef ENABLE_DRY_MODE
            !$omp do private(n, p, l, points, btot, q_c) &
            !$omp& private(is, js, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, dbuoyg, humg, tbuoyg, volg)
#else
            !$omp do private(n, p, l, points, btot) &
            !$omp& private(is, js, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, tbuoyg, volg)
#endif

            do n = 1, surf_parcels%local_num

#ifndef ENABLE_DRY_MODE
                ! liquid water content
                q_c = max(zero, surf_parcels%humidity(n) -  qq)

                ! total buoyancy (including effects of latent heating)
                btot = surf_parcels%buoyancy(n) + glat * q_c
#else
                btot = surf_parcels%buoyancy(n)
#endif


#ifdef ENABLE_BUOYANCY_PERTURBATION_MODE
                ! remove basic state N^2 * z
                btot = btot - bfsq * zpos
#endif

#ifndef ENABLE_P2G_1POINT
                points = surf_parcels%get_points(n)
#else
                points(:, 1) = surf_parcels%position(:, n)
#endif

!                 call get_index(surf_parcels%position(:, n), i, j)
!                 i = mod(i + nx, nx)
!                 j = mod(j + ny, ny)
!                 nparg(j, i) = nparg(j, i) + 1
!                 if (surf_parcels%area(n) <= amin) then
!                     nsparg(j, i) = nsparg(j, i) + 1
!                 endif

                btot = surf_parcels%buoyancy(n)

                ! we have 2 points per ellipse
                do p = 1, n_points_p2g

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    weights = point_weight_p2g * weights * surf_parcels%area(n)

                    do l = 1, 3
                        vortg(iz, js:js+1, is:is+1, l) = vortg(iz, js:js+1, is:is+1, l) &
                                                       + weights * surf_parcels%vorticity(l, n)
                    enddo
#ifndef ENABLE_DRY_MODE
                    dbuoyg(iz, js:js+1, is:is+1) = dbuoyg(iz, js:js+1, is:is+1) &
                                                 + weights * surf_parcels%buoyancy(n)
                    humg(iz, js:js+1, is:is+1) = humg(iz, js:js+1, is:is+1) &
                                               + weights * surf_parcels%humidity(n)
#endif
                    tbuoyg(iz, js:js+1, is:is+1) = tbuoyg(iz, js:js+1, is:is+1) &
                                                 + weights * btot

                    volg(iz, js:js+1, is:is+1) = volg(iz, js:js+1, is:is+1) &
                                               + weights
                enddo
            enddo
            !$omp end do
            !$omp end parallel

!             ! sanity check
!             if (sum(nparg) /= n_par) then
!                 print *, "par2grid: Wrong total number of parcels!"
!                 stop
!             endif

        end subroutine m_par2grid

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine surf_grid2par(add)
            logical, optional, intent(in) :: add

            call m_grid2par(0,  bot_parcels, add)
            call m_grid2par(nz, top_parcels, add)

        end subroutine surf_grid2par

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Interpolate the gridded quantities to the parcels
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        ! @pre The parcel must be assigned to the correct MPI process and the halo of fields must be
        !      filled correctly.
        subroutine m_grid2par(iz, surf_parcels, add)
            integer,               intent(in)    :: iz
            type(ellipse_pc_type), intent(inout) :: surf_parcels
            logical, optional, intent(in)        :: add
            double precision                     :: points(2, 2), weights(0:1, 0:1), dvdx(0:1, 0:1)
            integer                              :: n, p, l, is, js

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, surf_parcels%local_num
                        surf_parcels%delta_pos(:, n) = zero
                        surf_parcels%delta_vor(:, n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, surf_parcels%local_num
                    surf_parcels%delta_pos(:, n) = zero
                    surf_parcels%delta_vor(:, n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, is, js, weights, dvdx)
            do n = 1, surf_parcels%local_num

                surf_parcels%strain(:, n) = zero

#ifndef ENABLE_G2P_1POINT
                points = surf_parcels%get_points(n)
#else
                points(:, 1) = parcels%position(:, n)
#endif

                ! we have 2 points per ellipse
                do p = 1, n_points_g2p

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    do l = 1, 2
                        surf_parcels%delta_pos(:, n) = surf_parcels%delta_pos(:, n)                 &
                                                     + point_weight_g2p                             &
                                                     * sum(weights * velog(iz, js:js+1, is:is+1, l))
                    enddo

                    ! du/dx
                    surf_parcels%strain(1, n) = surf_parcels%strain(1, n)                       &
                                              + point_weight_g2p                                &
                                              * sum(weights * velgradg(iz, js:js+1, is:is+1, 1))

                    ! du/dy
                    surf_parcels%strain(2, n) = surf_parcels%strain(2, n)                       &
                                              + point_weight_g2p                                &
                                              * sum(weights * velgradg(iz, js:js+1, is:is+1, 2))

                    ! dv/dx = \zeta + du/dy
                    dvdx = vortg(iz, js:js+1, is:is+1, 3) + velgradg(iz, js:js+1, is:is+1, 2)

                    surf_parcels%strain(3, n) = surf_parcels%strain(3, n)                       &
                                              + point_weight_g2p * sum(weights * dvdx)

                    ! dv/dy
                    surf_parcels%strain(4, n) = surf_parcels%strain(4, n)                       &
                                              + point_weight_g2p                                &
                                              * sum(weights * velgradg(iz, js:js+1, is:is+1, 3))

                    do l = 1, 3
                        surf_parcels%delta_vor(:, n) = surf_parcels%delta_vor(:, n)                 &
                                                     + point_weight_g2p                             &
                                                     + sum(weights * vtend(iz, js:js+1, is:is+1, l))
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine m_grid2par

end module parcel_ellipse_interpl
