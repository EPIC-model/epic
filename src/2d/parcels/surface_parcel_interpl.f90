! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation for surface parcels.
! =============================================================================
module surface_parcel_interpl
    use constants, only : zero, one, two
    use parameters, only : nx, nz, lmin
    use options, only : parcel
    use surface_parcel_container, only : top_parcels, n_top_parcels &
                                       , bot_parcels, n_bot_parcels
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

    private :: is, weights, _len2grid, _surf_par2grid

    contains

        subroutine len2grid

            call _len2grid(nz, n_top_parcels, top_parcels)
            call _len2grid(0,  n_bot_parcels, bot_parcels)

        end subroutine len2grid

        ! Interpolate the parcel length to the grid
        subroutine _len2grid(iz, n_par, spar)
            integer,                        intent(in)    :: iz
            integer,                        intent(in)    :: n_par
            type(surface_parcel_container), intent(inout) :: spar
            double precision                              :: points(2)
            integer                                       :: n, p, l

            volg(iz, :) = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pl, is, js, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_par
                points(1) = spar%position(n) + f14 * spar%length(n)
                points(2) = spar%position(n) - f14 * spar%length(n)

                ! we have 2 points per line
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_surface_periodic_bc(points(p))

                    ! get interpolation weights and mesh indices
                    call linear(points(p), is, weights)

                    do l = 1, ngp
                        volg(iz, is(l)) = volg(iz, is(l)) &
                                        + f12 * weights(l) * spar%length(n)
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine _len2grid


        subroutine surf_par2grid

            call _surf_par2grid(nz, n_top_parcels, top_parcels)
            call _surf_par2grid(0,  n_bot_parcels, bot_parcels)

        end subroutine surf_par2grid

        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        !   - length
        ! It also updates the scalar fields:
        subroutine _surf_par2grid(iz, n_par, spar)
            integer,                        intent(in)    :: iz
            integer,                        intent(in)    :: n_par
            type(surface_parcel_container), intent(inout) :: spar
            double precision                              :: points(2)
            integer                                       :: n, p, l, i
            double precision                              :: weight, btot
#ifndef ENABLE_DRY_MODE
            double precision                              :: q_c
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
            !$omp do private(n, p, l, i, points, weight, btot, q_c, is, weights) &
            !$omp& reduction(+:vortg, dbuoyg, humg, tbuoyg, volg)
#else
            !$omp do private(n, p, l, i, points, weight, btot, is, weights) &
            !$omp& reduction(+:vortg, tbuoyg, volg)
#endif
            do n = 1, n_par

#ifndef ENABLE_DRY_MODE
                ! liquid water content
                q_c = parcels%humidity(n) &
                    - q_0 * dexp(lambda_c * (lower(2) - parcels%position(2, n)))
                q_c = max(zero, q_c)

                ! total buoyancy (including effects of latent heating)
                btot = parcels%buoyancy(n) + glat * q_c
#else
                btot = parcels%buoyancy(n)
#endif

                points(1) = spar%position(n) + f14 * spar%length(n)
                points(2) = spar%position(n) - f14 * spar%length(n)

                ! we have 2 points per line
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_surface_periodic_bc(points(p))

                    ! get interpolation weights and mesh indices
                    call linear(points(p), is, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per line
                    do l = 1, ngp

                        weight = f12 * weights(l) * parcels%length(n)

                        vortg(iz, is(l)) = vortg(iz, is(l)) &
                                         + weight * parcels%vorticity(n)

#ifndef ENABLE_DRY_MODE
                        dbuoyg(iz, is(l)) = dbuoyg(iz, is(l)) &
                                          + weight * parcels%buoyancy(n)
                        humg(iz, is(l)) = humg(iz, is(l)) &
                                        + weight * parcels%humidity(n)
#endif
                        tbuoyg(iz, is(l)) = tbuoyg(iz, is(l)) + weight * btot
                        volg(iz, is(l)) = volg(iz, is(l)) + weight
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine _surf_par2grid


        ! Interpolate the gridded quantities to the parcels
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vor is the parcel vorticity
        ! @param[inout] vgrad is the parcel strain
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine surf_grid2par(iz, n_par, spar, vel, vor, vgrad, add)
            integer,                        intent(in)    :: iz
            integer,                        intent(in)    :: n_par
            type(surface_parcel_container), intent(inout) :: spar
            double precision,               intent(inout) :: vel, vor, vgrad(:)
            logical, optional,              intent(in)    :: add
            double precision                              :: points(2), weight
            integer                                       :: n, p, l

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_par
                        vel(n) = zero
                        vor(n)    = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_par
                    vel(n) = zero
                    vor(n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, weight, is, weights)
            do n = 1, n_par

                vgrad(n) = zero

                points(1) = spar%position(n) + f14 * spar%length(n)
                points(2) = spar%position(n) - f14 * spar%length(n)

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
                        vel(n) = vel(n) &
                               + weight * velog(iz, is(l), :)

                        vgrad(n) = vgrad(n) &
                                 + weight * velgradg(iz, is(l), :)

                        vor(n) = vor(n) + weight * vtend(iz, is(l))
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine surf_grid2par


        ! Interpolate the gridded quantities to the parcels without resetting
        ! their values to zero before doing grid2par.
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vor is the parcel vorticity
        ! @param[inout] vgrad is the parcel strain
        subroutine surf_grid2par_add(iz, n_par, spar, vel, vor, vgrad, add)
            integer,                        intent(in)    :: iz
            integer,                        intent(in)    :: n_par
            type(surface_parcel_container), intent(inout) :: spar
            double precision,               intent(inout) :: vel, vor, vgrad(:)

            call surf_grid2par(iz, n_par, spar, vel, vor, vgrad, add=.true.)

        end subroutine surf_grid2par_add


        ! Linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] ww interpolation weights
        subroutine linear(pos, ii, ww)
            double precision, intent(in)  :: pos
            integer,          intent(out) :: ii(2)
            double precision, intent(out) :: ww(2)
            double precision              :: x

            ! (i)
            ii(1) = floor((pos - lower(1)) * dxi(1))

            x = lower(1) + dble(i) * dx(1)

            ww(1) = one - abs(pos - x) * dxi(1)

            ! (i+1)
            ii(2) = ii(1) + 1
            x = lower(1) + dble(ii(2)) * dx(1)
            ww(2) = one - abs(pos - x) * dxi(1)

            ! account for x periodicity
            call periodic_index_shift(ii)

        end subroutine linear

end module surface_parcel_interpl
