! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : zero, one, two, f14
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipsoid
    use fields
    use phys_constants, only : h_0
    use phys_parameters, only : glat, lam_c
    use omp_lib
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = 8

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is(ngp), js(ngp), ks(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

    integer :: par2grid_timer, &
#ifndef NDBEBUG
               sym_vol2grid_timer, &
#endif
               grid2par_timer

    private :: is, js, ks, weights

    contains

        ! Interpolate the parcel volume to the grid
        subroutine vol2grid
            double precision  :: points(3, 4)
            integer           :: n, p, l
            double precision  :: pvol

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, is, js, ks, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n))


                ! we have 4 points per ellipsoid
                do p = 1, 4

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(:, p), is, js, ks, weights)

                    do l = 1, ngp
                        volg(ks(l), js(l), is(l)) = volg(ks(l), js(l), is(l)) &
                                                  + f14 * weights(l) * pvol
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

#ifndef NDEBUG
        ! Interpolate the parcel volume to the grid to check symmetry
        subroutine vol2grid_symmetry_error
            double precision :: points(3, 4), V, B(5), pos(3)
            integer          :: n, p, l, m
            double precision :: pvol

            call start_timer(sym_vol2grid_timer)

            sym_volg = zero

            do m = -1, 1, 2
                !$omp parallel default(shared)
                !$omp do private(n, p, l, points, pos, pvol, V, B, is, js, ks, weights) &
                !$omp& reduction(+: sym_volg)
                do n = 1, n_parcels

                    pos = parcels%position(:, n)
                    pvol = parcels%volume(n)
                    pos(1) = dble(m) * pos(1)
                    V = dble(m) * pvol
                    B = parcels%B(:, n)

                    B(2) = dble(m) * B(2)

                    points = get_ellipsoid_points(pos, V, B)

                    ! we have 4 points per ellipsoid
                    do p = 1, 4

                        ! ensure point is within the domain
                        call apply_periodic_bc(points(:, p))

                        ! get interpolation weights and mesh indices
                        call trilinear(points(:, p), is, js, ks, weights)

                        do l = 1, ngp
                            sym_volg(ks(l), js(l), is(l)) = sym_volg(ks(l), js(l), is(l)) &
                                                          + dble(m) * f14 * weights(l) * pvol
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
            double precision :: points(3, 4)
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
                ! liquid water content
                h_c = parcels%humidity(n) &
                    - h_0 * dexp(lam_c * (lower(3) - parcels%position(3, n)))
                h_c = max(zero, h_c)

                ! total buoyancy (including effects of latent heating)
                btot = parcels%buoyancy(n) + glat * h_c
#else
                btot = parcels%buoyancy(n)
#endif
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n))

                call get_index(parcels%position(:, n), i, j, k)
                i = mod(i + nx, nx)
                j = mod(j + ny, ny)
                nparg(k, j, i) = nparg(k, j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(k, j, i) = nsparg(k, j, i) + 1
                endif

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is a quarter due to 4 points per ellipsoid
                    do l = 1, ngp

                        weight = f14 * weights(l) * pvol

                        vortg(ks(l), js(l), is(l), :) = vortg(ks(l), js(l), is(l), :) &
                                                      + weight * parcels%vorticity(:, n)

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
            dbuoyg(0,  :, :) = two * dbuoyg(0,  :, :)
            dbuoyg(nz, :, :) = two * dbuoyg(nz, :, :)
            dbuoyg(1,    :, :) = dbuoyg(1,    :, :) + dbuoyg(-1,   :, :)
            dbuoyg(nz-1, :, :) = dbuoyg(nz-1, :, :) + dbuoyg(nz+1, :, :)
#endif
            tbuoyg(0,  :, :) = two * tbuoyg(0,  :, :)
            tbuoyg(nz, :, :) = two * tbuoyg(nz, :, :)
            tbuoyg(1,    :, :) = tbuoyg(1,    :, :) + tbuoyg(-1,   :, :)
            tbuoyg(nz-1, :, :) = tbuoyg(nz-1, :, :) + tbuoyg(nz+1, :, :)

            ! exclude halo cells to avoid division by zero
            do p = 1, 3
                vortg(0:nz, :, :, p) = vortg(0:nz, :, :, p) / volg(0:nz, :, :)
            enddo

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
            double precision                    :: points(3, 4), weight
            integer                             :: n, p, c, l

            call start_timer(grid2par_timer)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_parcels
                        vel(:, n) = zero
                        vor(:, n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    vel(:, n) = zero
                    vor(:, n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, c, points, weight, is, js, ks, weights)
            do n = 1, n_parcels

                vgrad(:, n) = zero

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              parcels%volume(n),      &
                                              parcels%B(:, n))

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        weight = f14 * weights(l)

                        ! loop over field components
                        do c = 1, 3
                            ! the weight is a quarter due to 4 points per ellipsoid
                            vel(c, n) = vel(c, n) &
                                      + weight * velog(ks(l), js(l), is(l), c)
                        enddo

                        do c = 1, 5
                            vgrad(c, n) = vgrad(c, n) &
                                        + weight * velgradg(ks(l), js(l), is(l), c)
                        enddo

                        vor(:, n) = vor(:, n) + weight * vtend(ks(l), js(l), is(l), :)
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


        ! Tri-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii zonal grid points for interoplation
        ! @param[out] jj meridional grid points for interpolation
        ! @param[out] kk vertical grid points for interpolation
        ! @param[out] ww interpolation weights
        subroutine trilinear(pos, ii, jj, kk, ww)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: ii(ngp), jj(ngp), kk(ngp)
            double precision, intent(out) :: ww(ngp)
            double precision              :: xyz(3)
            integer                       :: i, l

            ! (i, j, k)
            call get_index(pos, ii(1), jj(1), kk(1))

            ! (i, j, k+1)
            ii(5) = ii(1)
            jj(5) = jj(1)
            kk(5) = kk(1) + 1

            do i = 0, 1
                l = 4 * i
                call get_position(ii(1+l), jj(1+l), kk(1+l), xyz)
                ww(1+l) = product(one - abs(pos - xyz) * dxi)

                ! (i+1, j, k) and (i+1, j, k+1)
                ii(2+l) = ii(1+l) + 1
                jj(2+l) = jj(1+l)
                kk(2+l) = kk(1+l)
                call get_position(ii(2+l), jj(2+l), kk(2+l), xyz)
                ww(2+l) = product(one - abs(pos - xyz) * dxi)

                ! (i, j+1, k) and (i, j+1, k+1)
                ii(3+l) = ii(1+l)
                jj(3+l) = jj(1+l) + 1
                kk(3+l) = kk(1+l)
                call get_position(ii(3+l), jj(3+l), kk(3+l), xyz)
                ww(3+l) = product(one - abs(pos - xyz) * dxi)

                ! (i+1, j+1, k) an (i+1, j+1, k+1)
                ii(4+l) = ii(2+l)
                jj(4+l) = jj(3+l)
                kk(4+l) = kk(1+l)
                call get_position(ii(4+l), jj(4+l), kk(4+l), xyz)
                ww(4+l) = product(one - abs(pos - xyz) * dxi)
            enddo

            ! account for x and y periodicity
            call periodic_index_shift(ii, jj)

        end subroutine trilinear

end module parcel_interpl
