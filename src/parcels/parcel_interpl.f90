! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : zero, one, two, ndim
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, ny, nz, vmin
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipse
    use fields
    use phys_constants, only : h_0
    use phys_parameters, only : glat, lam_c
    use omp_lib
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = ndim ** 2

    ! number of ellipse / ellispoid support points
    integer, parameter :: nsup = 2 ** (ndim - 1)

    ! interpolation indices
    integer :: is(ngp), js(ngp), ks(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

    integer :: par2grid_timer, &
#ifndef NDBEBUG
               sym_vol2grid_timer, &
#endif
               grid2par_timer


    private :: is, js, ks, weights, nsup

    contains

        ! Interpolate the parcel volume to the grid
        subroutine vol2grid
            double precision  :: points(nsup, ndim)
            integer           :: n, p, l
            double precision  :: pvol

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, is, js, ks, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(n, :), &
                                            pvol, parcels%B(n, :))


                ! iterate over support points (2D --> 2, 3D --> 4)
                do p = 1, nsup

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call linterp(points(p, :), weights, is, js, ks)

                    do l = 1, ngp
                        volg(ks(l), js(l), is(l)) = volg(ks(l), js(l), is(l)) &
                                                  + f12 * weights(l) * pvol
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ! apply free slip boundary condition
            volg(0,  :, :) = two * volg(0,  :, :)
            volg(nz, :, :) = two * volg(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :, :) = volg(1,    :, :) + volg(-1,   :, :)
            volg(nz-1, :, :) = volg(nz-1, :, :) + volg(nz+1, :, :)

        end subroutine vol2grid

        !FIXME
#ifndef NDEBUG
        ! Interpolate the parcel volume to the grid to check symmetry
        subroutine vol2grid_symmetry_error
            double precision :: points(nsup, ndim), V, B(2), pos(2)
            integer          :: n, p, l, m
            double precision :: pvol

            call start_timer(sym_vol2grid_timer)

            sym_volg = zero

            do m = -1, 1, 2
                !$omp parallel default(shared)
                !$omp do private(n, p, l, points, pos, pvol, V, B, is, js, ks, weights) &
                !$omp& reduction(+: sym_volg)
                do n = 1, n_parcels

                    pos = parcels%position(n, :)
                    pvol = parcels%volume(n)
                    pos(1) = dble(m) * pos(1)
                    V = dble(m) * pvol
                    B = parcels%B(n, :)

                    B(2) = dble(m) * B(2)

                    points = get_ellipse_points(pos, V, B)

                    ! we have 2 points per ellipse
                    do p = 1, 2

                        ! ensure point is within the domain
                        call apply_periodic_bc(points(p, :))

                        ! get interpolation weights and mesh indices
                        call linterp(points(p, :), weights, is, js, ks)

                        do l = 1, ngp
                            sym_volg(ks(l), js(l), is(l)) = sym_volg(ks(l), js(l), is(l)) &
                                                          + dble(m) * f12 * weights(l) * pvol
                        enddo
                    enddo
                enddo
                !$omp end do
                !$omp end parallel
            enddo
            call stop_timer(sym_vol2grid_timer)
        end subroutine vol2grid_symmetry_error
#endif

        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        !   - volume
        ! It also updates the scalar fields:
        !   - nparg, that is the number of parcels per grid cell
        !   - nsparg, that is the number of small parcels per grid cell
        subroutine par2grid
            double precision :: points(nsup, ndim)
            integer          :: n, p, l, i, j, k
            double precision :: pvol, weight, btot
#ifndef ENABLE_DRY_MODE
            double precision :: h_c
#endif

            call start_timer(par2grid_timer)

            vortg = zero
            volg = zero
            nparg = zero
            nsparg = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg = zero
#endif
            tbuoyg = zero
            !$omp parallel default(shared)
#ifndef ENABLE_DRY_MODE
            !$omp do private(n, p, l, i, j, k, points, pvol, weight, btot, h_c, is, js, ks, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, dbuoyg, tbuoyg, volg)
#else
            !$omp do private(n, p, l, i, j, k, points, pvol, weight, btot, is, js, ks, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, tbuoyg, volg)
#endif
            do n = 1, n_parcels
                pvol = parcels%volume(n)

#ifndef ENABLE_DRY_MODE
                !FIXME
                ! liquid water content
                h_c = parcels%humidity(n) &
                    - h_0 * dexp(lam_c * (lower(2) - parcels%position(n, 2)))
                h_c = max(zero, h_c)

                ! total buoyancy (including effects of latent heating)
                btot = parcels%buoyancy(n) + glat * h_c
#else
                btot = parcels%buoyancy(n)
#endif
                points = get_ellipse_points(parcels%position(n, :), &
                                            pvol, parcels%B(n, :))

                call get_index(parcels%position(n, :), i, j)
                i = mod(i + nx, nx)
                nparg(k, j, i) = nparg(k, j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(k, j, i) = nsparg(k, j, i) + 1
                endif

                ! iterate over support points (2D --> 2, 3D --> 4)
                do p = 1, nsup

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call linterp(points(p, :), weights, is, js, ks)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    do l = 1, ngp

                        weight = f12 * weights(l) * pvol

                        vortg(ks(l), js(l), is(l), :) = vortg(ks(l), js(l), is(l), :) &
                                                      + weight * parcels%vorticity(n, :)

#ifndef ENABLE_DRY_MODE
                        dbuoyg(ks(l), js(l), is(l)) = dbuoyg(ks(l), js(l), is(l)) &
                                                    + weight * parcels%buoyancy(n)
#endif
                        tbuoyg(ks(l), js(l), is(l)) = tbuoyg(ks(l), js(l), is(l)) &
                                                    + weight * btot
                        volg(ks(l), js(l), is(l)) = volg(ks(l), js(l), is(l)) &
                                                  + weight
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ! apply free slip boundary condition
            volg(0,  :, :) = two * volg(0,  :, :)
            volg(nz, :, :) = two * volg(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :, :) = volg(1,    :, :) + volg(-1,   :, :)
            volg(nz-1, :, :) = volg(nz-1, :, :) + volg(nz+1, :, :)

            vortg(0,  :, :, :) = two * vortg(0,  :, :, :)
            vortg(nz, :, :, :) = two * vortg(nz, :, :, :)
            vortg(1,    :, :, :) = vortg(1,    :, :, :) + vortg(-1,   :, :, :)
            vortg(nz-1, :, :, :) = vortg(nz-1, :, :, :) + vortg(nz+1, :, :, :)

#ifndef ENABLE_DRY_MODE
            dbuoyg(0,  :, :) = two * dbuoyg(0,  :)
            dbuoyg(nz, :, :) = two * dbuoyg(nz, :)
            dbuoyg(1,    :, :) = dbuoyg(1,    :, :) + dbuoyg(-1,   :, :)
            dbuoyg(nz-1, :, :) = dbuoyg(nz-1, :, :) + dbuoyg(nz+1, :, :)
#endif
            tbuoyg(0,  :, :) = two * tbuoyg(0,  :, :)
            tbuoyg(nz, :, :) = two * tbuoyg(nz, :, :)
            tbuoyg(1,    :, :) = tbuoyg(1,    :, :) + tbuoyg(-1,   :, :)
            tbuoyg(nz-1, :, :) = tbuoyg(nz-1, :, :) + tbuoyg(nz+1, :, :)

            ! exclude halo cells to avoid division by zero
            do i = 1, vdim
                vortg(0:nz, :, :, i) = vortg(0:nz, :, :, i) / volg(0:nz, :, :)
            enddo

            ! extrapolate to halo grid points (since halo grid points
            ! are used to get u_z = w_x - zeta)
            vortg(-1,   :, :, :) = two * vortg(0,  :, :, :) - vortg(1,    :, :, :)
            vortg(nz+1, :, :, :) = two * vortg(nz, :, :, :) - vortg(nz-1, :, :, :)

#ifndef ENABLE_DRY_MODE
            dbuoyg(0:nz, :, :) = dbuoyg(0:nz, :, :) / volg(0:nz, :, :)
#endif
            tbuoyg(0:nz, :, :) = tbuoyg(0:nz, :, :) / volg(0:nz, :, :)

            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            tbuoyg(-1,   :, :) = two * tbuoyg(0,  :, :) - tbuoyg(1, :, :)
            tbuoyg(nz+1, :, :) = two * tbuoyg(nz, :, :) - tbuoyg(nz-1, :, :)

            ! sum halo contribution into internal cells
            ! (be aware that halo cell contribution at upper boundary
            ! are added to cell nz)
            nparg(0,    :, :) = nparg(0,    :, :) + nparg(-1, :, :)
            nparg(nz-1, :, :) = nparg(nz-1, :, :) + nparg(nz, :, :)

            nsparg(0,    :, :) = nsparg(0,    :, :) + nsparg(-1, :, :)
            nsparg(nz-1, :, :) = nsparg(nz-1, :, :) + nsparg(nz, :, :)

            ! sanity check
            if (sum(nparg(0:nz-1, :, :)) /= n_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

            call stop_timer(par2grid_timer)

        end subroutine par2grid


        ! Interpolate the gridded quantities to the parcels
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vor is the parcel vorticity
        ! @param[inout] vgrad is the parcel strain
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine grid2par(vel, vor, vgrad, add)
            double precision,     intent(inout) :: vel(:, :), vor(:, :), vgrad(:, :)
            logical, optional, intent(in)       :: add
            double precision                    :: points(nsup, ndim), weight
            integer                             :: n, p, c, l

            call start_timer(grid2par_timer)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_parcels
                        vel(n, :) = zero
                        vor(n, :) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    vel(n, :) = zero
                    vor(n, :) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, c, points, weight, is, js, ks, weights)
            do n = 1, n_parcels

                vgrad(n, :) = zero

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n),      &
                                            parcels%B(n, :))

                ! iterate over support points (2D --> 2, 3D --> 4)
                do p = 1, nsup

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call linterp(points(p, :), weights, is, js, ks)

                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        weight = f12 * weights(l)

                        ! loop over field components
                        do c = 1, ndim
                            ! the weight is halved due to 2 points per ellipse
                            vel(n, c) = vel(n, c) &
                                      + weight * velog(ks(l), js(l), is(l), c)
                        enddo

                        do c = 1, 4 !FIXME
                            vgrad(n, c) = vgrad(n, c) &
                                        + weight * velgradg(ks(l), js(l), is(l), c)
                        enddo

                        do c = 1, vdim
                            vor(n, c) = vor(n, c) + weight * vtend(ks(l), js(l), is(l))
                        enddo
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(grid2par_timer)

        end subroutine grid2par


        ! Interpolate the gridded quantities to the parcels without resetting
        ! their values to zero before doing grid2par.
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vor is the parcel vorticity
        ! @param[inout] vgrad is the parcel strain
        subroutine grid2par_add(vel, vor, vgrad)
            double precision,       intent(inout) :: vel(:, :), vor(:, :), vgrad(:, :)

            call grid2par(vel, vor, vgrad, add=.true.)

        end subroutine grid2par_add


        ! Linear interpolation: bi-linear in 2D, tri-linear in 3D
        ! @param[in] pos position of the parcel
        ! @param[out] ww interpolation weights
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj longitudinal (2D: vertical) grid points for interpolation
        ! @param[out] kk vertical (3D only) grid points for interpolation
        subroutine linterp(pos, ww, ii, jj, kk)
            double precision,  intent(in)  :: pos(ndim)
            double precision,  intent(out) :: ww(ngp)
            integer,           intent(out) :: ii(ngp), jj(ngp), kk(ngp)

#ifdef ENABLE_3D
            call trilinear(pos, ww, ii, jj, kk)
#else
            call bilinear(pos, ww, ii, jj)
#endif
        end subroutine linterp

#ifdef ENABLE_3D
        ! Tri-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ww interpolation weights
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj vertical grid points for interpolation
        subroutine trilinear(pos, ww, ii, jj, kk)
            double precision, intent(in)  :: pos(ndim)
            double precision, intent(out) :: ww(ngp)
            integer,          intent(out) :: ii(ngp), jj(ngp), kk(ngp)
            double precision              :: xyz(ndim)
            integer                       :: k, j

            ! (i, j, k)
            call get_index(pos, ii(1), jj(1), kk(1))

            kk(5) = kk(1) + 1

            do k = 0, 1
                j = 4 * k
                ! (i, j, k) and (i, j, k+1)
                call get_position((/ii(1 + j), jj(1 + j), kk(1 + j)/), xyz)
                ww(1 + j) = product(one - abs(pos - xyz) * dxi)

                ! (i+1, j, k) and (i+1, j, k+1)
                ii(2 + j) = ii(1 + j) + 1
                jj(2 + j) = jj(1 + j)
                kk(2 + j) = kk(1 + j)
                call get_position((/ii(2 + j), jj(2 + j), kk(2 + j)/), xyz)
                ww(2 + j) = product(one - abs(pos - xyz) * dxi)

                ! (i, j+1, k) and (i, j+1, k+1)
                ii(3 + j) = ii(1 + j)
                jj(3 + j) = jj(1 + j) + 1
                kk(3 + j) = kk(1 + j)
                call get_position((/ii(3 + j), jj(3 + j), kk(3 + j)/), xyz)
                ww(3 + j) = product(one - abs(pos - xyz) * dxi)

                ! (i+1, j+1, k) and (i+1, j+1, k+1)
                ii(4 + j) = ii(2 + j)
                jj(4 + j) = jj(3 + j)
                kk(4 + j) = kk(1 + j)
                call get_position((/ii(4 + j), jj(4 + j), kk(4 + j)/), xyz)
                ww(4 + j) = product(one - abs(pos - xyz) * dxi)
            enddo

            ! account for x periodicity
            call periodic_index_shift(ii, nx)

            ! account for y periodicity
            call periodic_index_shift(jj, ny)

        end subroutine trilinear
#else
        ! Bi-linear
        ! @param[in] pos position of the parcel
        ! @param[out] ww interpolation weights
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj vertical grid points for interpolation
        subroutine bilinear(pos, ww, ii, jj)
            double precision, intent(in)  :: pos(ndim)
            double precision, intent(out) :: ww(ngp)
            integer,          intent(out) :: ii(ngp), jj(ngp)
            double precision              :: xy(ndim)

            ! (i, j)
            call get_index(pos, ii(1), jj(1))
            call get_position((/ii(1), jj(1)/), xy)
            ww(1) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j)
            ii(2) = ii(1) + 1
            jj(2) = jj(1)
            call get_position((/ii(2), jj(2)/), xy)
            ww(2) = product(one - abs(pos - xy) * dxi)

            ! (i, j+1)
            ii(3) = ii(1)
            jj(3) = jj(1) + 1
            call get_position((/ii(3), jj(3)/), xy)
            ww(3) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j+1)
            ii(4) = ii(2)
            jj(4) = jj(3)
            call get_position((/ii(4), jj(4)/), xy)
            ww(4) = product(one - abs(pos - xy) * dxi)

            ! account for x periodicity
            call periodic_index_shift(ii, nx)

        end subroutine bilinear
#endif

end module parcel_interpl
