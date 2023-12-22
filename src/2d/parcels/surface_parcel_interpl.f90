! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation for surface parcels.
! =============================================================================
module surface_parcel_interpl
    use constants, only : zero, one, f12, f14, f34
    use parameters, only : nx, nz, lmin
    use options, only : parcel
    use surface_parcel_container, only : top_parcels, n_top_parcels     &
                                       , bot_parcels, n_bot_parcels     &
                                       , surface_parcel_container_type  &
                                       , get_surface_parcel_length
    use surface_parcel_bc, only : apply_surface_periodic_bc
    use fields
    use physics, only : glat, lambda_c, q_0
    use omp_lib
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = 2

    ! interpolation indices
    integer :: is(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

    private :: is, weights, len2grid_, surf_par2grid_, surf_grid2par_, get_line_points

    contains

        function get_line_points(n, spar) result(points)
            integer,                             intent(in) :: n
            type(surface_parcel_container_type), intent(in) :: spar
            double precision                                :: length
            double precision                                :: points(2)

            length = get_surface_parcel_length(n, spar)

            points(1) = spar%position(n) + f14 * length
            points(2) = spar%position(n) + f34 * length

        end function get_line_points

        subroutine len2grid

            call len2grid_(nz, n_top_parcels, top_parcels)
            call len2grid_(0,  n_bot_parcels, bot_parcels)

        end subroutine len2grid

        ! Interpolate the parcel length to the grid
        subroutine len2grid_(iz, n_par, spar)
            integer,                             intent(in) :: iz
            integer,                             intent(in) :: n_par
            type(surface_parcel_container_type), intent(in) :: spar
            double precision                                :: points(2), length
            integer                                         :: n, p, l

            volg(iz, :) = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, is, weights, length) &
            !$omp& reduction(+: volg)
            do n = 1, n_par

                length = get_surface_parcel_length(n, spar)

                points = get_line_points(n, spar)


                ! we have 2 points per line
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_surface_periodic_bc(points(p))

                    ! get interpolation weights and mesh indices
                    call linear(points(p), is, weights)

                    do l = 1, ngp
                        volg(iz, is(l)) = volg(iz, is(l)) &
                                        + f12 * weights(l) * length
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine len2grid_


        subroutine surf_par2grid

            call surf_par2grid_(nz, n_top_parcels, top_parcels)
            call surf_par2grid_(0,  n_bot_parcels, bot_parcels)

        end subroutine surf_par2grid

        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        ! It also updates the scalar fields:
        subroutine surf_par2grid_(iz, n_par, spar)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            double precision                                   :: points(2)
            integer                                            :: n, p, l
            double precision                                   :: weight, btot, length
#ifndef ENABLE_DRY_MODE
            double precision                                   :: q_c
#endif
            vortg(iz, :) = zero
            volg(iz, :) = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg(iz, :) = zero
            humg(iz, :) = zero
#endif
            tbuoyg(iz, :) = zero
            !$omp parallel default(shared)
#ifndef ENABLE_DRY_MODE
            !$omp do private(n, p, l, points, weight, btot, q_c, is, weights, length) &
            !$omp& reduction(+:vortg, dbuoyg, humg, tbuoyg, volg)
#else
            !$omp do private(n, p, l, points, weight, btot, is, weights, length) &
            !$omp& reduction(+:vortg, tbuoyg, volg)
#endif
            do n = 1, n_par

#ifndef ENABLE_DRY_MODE
                ! liquid water content
                q_c = spar%humidity(n) &
                    - q_0 * dexp(- lambda_c * dble(iz) * dx(2))
                q_c = max(zero, q_c)

                ! total buoyancy (including effects of latent heating)
                btot = spar%buoyancy(n) + glat * q_c
#else
                btot = spar%buoyancy(n)
#endif

                length = get_surface_parcel_length(n, spar)

                points = get_line_points(n, spar)

                ! we have 2 points per line
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_surface_periodic_bc(points(p))

                    ! get interpolation weights and mesh indices
                    call linear(points(p), is, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per line
                    do l = 1, ngp

                        weight = f12 * weights(l) * length

                        vortg(iz, is(l)) = vortg(iz, is(l)) &
                                         + weight * spar%vorticity(n)

#ifndef ENABLE_DRY_MODE
                        dbuoyg(iz, is(l)) = dbuoyg(iz, is(l)) &
                                          + weight * spar%buoyancy(n)
                        humg(iz, is(l)) = humg(iz, is(l)) &
                                        + weight * spar%humidity(n)
#endif
                        tbuoyg(iz, is(l)) = tbuoyg(iz, is(l)) + weight * btot
                        volg(iz, is(l)) = volg(iz, is(l)) + weight
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine surf_par2grid_


        subroutine surf_grid2par(add)
            logical, optional, intent(in)    :: add

            call surf_grid2par_(0,  n_bot_parcels, bot_parcels, add)
            call surf_grid2par_(nz, n_top_parcels, top_parcels, add)

        end subroutine surf_grid2par


        ! Interpolate the gridded quantities to the parcels
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine surf_grid2par_(iz, n_par, spar, add)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            logical, optional,                   intent(in)    :: add
            double precision                                   :: points(2), weight
            integer                                            :: n, p, l

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_par
                        spar%delta_pos(n) = zero
                        spar%delta_vor(n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_par
                    spar%delta_pos(n) = zero
                    spar%delta_vor(n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, weight, is, weights)
            do n = 1, n_par

                points = get_line_points(n, spar)

                ! we have 2 points per line
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_surface_periodic_bc(points(p))

                    ! get interpolation weights and mesh indices
                    call linear(points(p), is, weights)

                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        weight = f12 * weights(l)

                        ! the weight is halved due to 2 points per line
                        spar%delta_pos(n) = spar%delta_pos(n) &
                                          + weight * velog(iz, is(l), 1)

                        spar%delta_vor(n) = spar%delta_vor(n) + weight * vtend(iz, is(l))
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine surf_grid2par_


        ! Linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] ww interpolation weights
        subroutine linear(pos, ii, ww)
            double precision, intent(in)  :: pos
            integer,          intent(out) :: ii(2)
            double precision, intent(out) :: ww(2)
            double precision              :: x, px

            ! (i)
            x = (pos - lower(1)) * dxi(1)
            ii(1) = floor(x)
            ii(2) = ii(1) + 1

            px = x - dble(ii(1))
            ww(1) = one - px
            ww(2) = px

            ! account for x periodicity
            call periodic_index_shift(ii)

        end subroutine linear

end module surface_parcel_interpl
