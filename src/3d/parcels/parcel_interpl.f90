! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : zero, one, two, f14, pi, f12, f32
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin, l_bndry_zeta_zero
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipsoid
    use fields
    use physics, only : glat, lambda_c, q_0
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
               grid2par_timer

    private :: is, js, ks, weights

    contains

        ! Interpolate the parcel volume to the grid
        subroutine vol2grid(l_reuse)
            logical, optional, intent(in) :: l_reuse
            double precision              :: points(3, 4)
            integer                       :: n, p, l
            double precision              :: pvol
            integer :: i, j, k, i_stored, j_stored, k_stored

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, is, js, ks, i, j, k, i_stored, j_stored, k_stored, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
               pvol = parcels%volume(n)

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n),  &
                                              n, l_reuse)


                ! we have 4 points per ellipsoid
                do p = 1, 4

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))
                    call get_index(points(:, p), i, j, k)

                    if (p == 1) then
                      call trilinear(points(:, p), is, js, ks, weights)
                    elseif ((i == i_stored) .and. (j == j_stored) .and. (k == k_stored)) then
                      ! if point is in same grid cell, just add to weights
                      call trilinear_weights_add(points(:, p), i, j, k, weights)
                    else
                      ! if point is in different grid cell, save previously stored
                      ! weights first
                      ! loop over grid points which are part of the interpolation
                      ! the weight is a quarter due to 4 points per ellipsoid
                      do l = 1, ngp
                          volg(ks(l), js(l), is(l)) = volg(ks(l), js(l), is(l)) &
                                                    + f14 * weights(l) * pvol
                      enddo
                      call trilinear(points(:, p), is, js, ks, weights)
                    endif
                    i_stored = i
                    j_stored = j
                    k_stored = k
                enddo
                ! save contibutions at end of points loop
                do l = 1, ngp
                    volg(ks(l), js(l), is(l)) = volg(ks(l), js(l), is(l)) &
                                              + f14 * weights(l) * pvol
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ! apply free slip boundary condition
            !$omp parallel workshare
            volg(0,  :, :) = two * volg(0,  :, :)
            volg(nz, :, :) = two * volg(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :, :) = volg(1,    :, :) + volg(-1,   :, :)
            volg(nz-1, :, :) = volg(nz-1, :, :) + volg(nz+1, :, :)

            ! Do volume extrapolation last, to compensate for effective location
            ! of parcels
            volg(0,  :, :) = f32 * volg(0,  :, :) - f12 * volg(1, :, :)
            volg(nz, :, :) = f32 * volg(nz, :, :) - f12 * volg(nz-1, :, :)            
            !$omp end parallel workshare
        end subroutine vol2grid


        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        !   - volume
        ! It also updates the scalar fields:
        !   - nparg, that is the number of parcels per grid cell
        !   - nsparg, that is the number of small parcels per grid cell
        subroutine par2grid(l_reuse)
            logical, optional :: l_reuse
            double precision  :: points(3, 4)
            integer           :: n, p, l, i, j, k
            double precision  :: pvol, weight, btot
#ifndef ENABLE_DRY_MODE
            double precision  :: q_c
#endif

            call start_timer(par2grid_timer)

            vortg = zero
            volg = zero
            nparg = zero
            nsparg = zero
#ifndef ENABLE_DRY_MODE
            dbuoyg = zero
            humg = zero
#endif
            tbuoyg = zero
            !$omp parallel default(shared)
#ifndef ENABLE_DRY_MODE
            !$omp do private(n, p, l, i, j, k, points, pvol, weight, btot, q_c) &
            !$omp& private( is, js, ks, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, dbuoyg, humg, tbuoyg, volg)
#else
            !$omp do private(n, p, l, i, j, k, points, pvol, weight, btot) &
            !$omp& private( is, js, ks, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, tbuoyg, volg)
#endif
            do n = 1, n_parcels
                pvol = parcels%volume(n)

#ifndef ENABLE_DRY_MODE
                ! liquid water content
                q_c = parcels%humidity(n) &
                    - q_0 * dexp(lambda_c * (lower(3) - parcels%position(3, n)))
                q_c = max(zero, q_c)

                ! total buoyancy (including effects of latent heating)
                btot = parcels%buoyancy(n) + glat * q_c
#else
                btot = parcels%buoyancy(n)
#endif
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n), n, l_reuse)

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
                    call get_index(points(:, p), i, j, k)

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
                        humg(ks(l), js(l), is(l)) = humg(ks(l), js(l), is(l)) &
                                                    + weight * parcels%humidity(n)
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

            !$omp parallel workshare
            ! apply free slip boundary condition
            volg(0,  :, :) = two * volg(0,  :, :)
            volg(nz, :, :) = two * volg(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :, :) = volg(1,    :, :) + volg(-1,   :, :)
            volg(nz-1, :, :) = volg(nz-1, :, :) + volg(nz+1, :, :)

            vortg(0,  :, :, :) = two * vortg(0,  :, :, :)
            vortg(nz, :, :, :) = two * vortg(nz, :, :, :)
            !$omp end parallel workshare

            !$omp parallel workshare
            vortg(1,    :, :, :) = vortg(1,    :, :, :) + vortg(-1,   :, :, :)
            vortg(nz-1, :, :, :) = vortg(nz-1, :, :, :) + vortg(nz+1, :, :, :)

#ifndef ENABLE_DRY_MODE
            dbuoyg(0,  :, :) = two * dbuoyg(0,  :, :)
            dbuoyg(nz, :, :) = two * dbuoyg(nz, :, :)
            dbuoyg(1,    :, :) = dbuoyg(1,    :, :) + dbuoyg(-1,   :, :)
            dbuoyg(nz-1, :, :) = dbuoyg(nz-1, :, :) + dbuoyg(nz+1, :, :)
            humg(0,  :, :) = two * humg(0,  :, :)
            humg(nz, :, :) = two * humg(nz, :, :)
            humg(1,    :, :) = humg(1,    :, :) + humg(-1,   :, :)
            humg(nz-1, :, :) = humg(nz-1, :, :) + humg(nz+1, :, :)
#endif
            tbuoyg(0,  :, :) = two * tbuoyg(0,  :, :)
            tbuoyg(nz, :, :) = two * tbuoyg(nz, :, :)
            tbuoyg(1,    :, :) = tbuoyg(1,    :, :) + tbuoyg(-1,   :, :)
            tbuoyg(nz-1, :, :) = tbuoyg(nz-1, :, :) + tbuoyg(nz+1, :, :)
            !$omp end parallel workshare

            ! exclude halo cells to avoid division by zero
            do p = 1, 3
                vortg(0:nz, :, :, p) = vortg(0:nz, :, :, p) / volg(0:nz, :, :)
            enddo

            ! At boundary extrapolation
            vortg(0,  :, :, :) = f32 * vortg(0,  :, :, :) - f12 * vortg(1, :, :, :)
            vortg(nz, :, :, :) = f32 * vortg(nz, :, :, :) - f12 * vortg(nz-1, :, :, :)            

            !-------------------------------------------------------
            ! Set zeta = 0 on the boundary if required:
            if (l_bndry_zeta_zero(1)) then
                vortg(0, :, :, I_Z) = zero
            endif

            if (l_bndry_zeta_zero(2)) then
                vortg(nz, :, :, I_Z) = zero
            endif

            !$omp parallel workshare
            vortg(-1,   :, :, :) = two * vortg(0,  :, :, :) - vortg(1, :, :, :)
            vortg(nz+1, :, :, :) = two * vortg(nz, :, :, :) - vortg(nz-1, :, :, :)

#ifndef ENABLE_DRY_MODE
            dbuoyg(0:nz, :, :) = dbuoyg(0:nz, :, :) / volg(0:nz, :, :)
            humg(0:nz, :, :) = humg(0:nz, :, :) / volg(0:nz, :, :)
#endif
            tbuoyg(0:nz, :, :) = tbuoyg(0:nz, :, :) / volg(0:nz, :, :)

            ! At boundary extrapolation
            humg(0,  :, :) = f32 * humg(0,  :, :) - f12 * humg(1, :, :)
            humg(nz, :, :) = f32 * humg(nz, :, :) - f12 * humg(nz-1, :, :)
            dbuoyg(0,  :, :) = f32 * dbuoyg(0,  :, :) - f12 * dbuoyg(1, :, :)
            dbuoyg(nz, :, :) = f32 * dbuoyg(nz, :, :) - f12 * dbuoyg(nz-1, :, :)
            tbuoyg(0,  :, :) = f32 * tbuoyg(0,  :, :) - f12 * tbuoyg(1, :, :)
            tbuoyg(nz, :, :) = f32 * tbuoyg(nz, :, :) - f12 * tbuoyg(nz-1, :, :)

            ! Do volume extrapolation only after all divisions by volg
            volg(0,  :, :) = f32 * volg(0,  :, :) - f12 * volg(1, :, :)
            volg(nz, :, :) = f32 * volg(nz, :, :) - f12 * volg(nz-1, :, :)            

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
            !$omp end parallel workshare

            ! sanity check
            if (sum(nparg(0:nz-1, :, :)) /= n_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

            call stop_timer(par2grid_timer)

        end subroutine par2grid


        ! Interpolate the gridded quantities to the parcels
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vortend is the parcel vorticity tendency
        ! @param[inout] vgrad is the parcel strain
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine grid2par(vel, vortend, vgrad, add)
          double precision,     intent(inout) :: vel(3, n_parcels), &
                                                 vortend(3, n_parcels), &
                                                 vgrad(5, n_parcels)
          logical, optional, intent(in)       :: add
          double precision :: points(3, 4)
            integer                             :: n, l, p
            !           double precision :: vsum

            call start_timer(grid2par_timer)

            ! clear old data efficiently
            if (present(add)) then
               if (add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_parcels
                        vel(:, n) = zero
                        vortend(:, n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    vel(:, n) = zero
                    vortend(:, n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, l, p, points, is, js, ks, weights)
            do n = 1, n_parcels

               vgrad(:, n) = zero

               points = get_ellipsoid_points(parcels%position(:, n), &
                                             parcels%volume(n), parcels%B(:, n), n)

               do p = 1, 4
                  ! ensure point is within the domain
                  call apply_periodic_bc(points(:, p))

                  ! get interpolation weights and mesh indices
                  call trilinear(points(:, p), is, js, ks, weights)

                  ! loop over grid points which are part of the interpolation
                  do l = 1, ngp
                     vel(:, n) = vel(:, n) + f14 * weights(l) * velog(ks(l), js(l), is(l), :)
                     vgrad(:, n) = vgrad(:, n) + f14 * weights(l) * velgradg(ks(l), js(l), is(l), :)
                     vortend(:, n) = vortend(:, n) + f14 * weights(l) * vtend(ks(l), js(l), is(l), :)
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
        ! @param[inout] vortend is the parcel vorticity tendency
        ! @param[inout] vgrad is the parcel strain
        subroutine grid2par_add(vel, vortend, vgrad)
            double precision, intent(inout) :: vel(:, :), vortend(:, :), vgrad(:, :)

            call grid2par(vel, vortend, vgrad, add=.true.)

        end subroutine grid2par_add


        ! Tri-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii zonal grid points for interoplation
        ! @param[out] jj meridional grid points for interpolation
        ! @param[out] kk vertical grid points for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine trilinear(pos, ii, jj, kk, ww)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: ii(ngp), jj(ngp), kk(ngp)
            double precision, intent(out) :: ww(ngp)
            double precision              :: xyz(3)

            ww = zero

            ! (i, j, k)
            xyz = (pos - lower) * dxi
            ii(1) = floor(xyz(1))
            jj(1) = floor(xyz(2))
            kk(1) = floor(xyz(3))

            call get_weights(xyz, ii(1), jj(1), kk(1), ww)

            ii(1) = mod(ii(1) + nx, nx)
            jj(1) = mod(jj(1) + ny, ny)

            ! (i+1, j, k)
            ii(2) = mod(ii(1) + 1 + nx, nx)
            jj(2) = jj(1)
            kk(2) = kk(1)

            ! (i, j+1, k)
            ii(3) = ii(1)
            jj(3) = mod(jj(1) + 1 + ny, ny)
            kk(3) = kk(1)

            ! (i+1, j+1, k)
            ii(4) = ii(2)
            jj(4) = jj(3)
            kk(4) = kk(1)

            ! (i, j, k+1)
            ii(5) = ii(1)
            jj(5) = jj(1)
            kk(5) = kk(1) + 1

            ! (i+1, j, k+1)
            ii(6) = ii(2)
            jj(6) = jj(1)
            kk(6) = kk(5)

            ! (i, j+1, k+1)
            ii(7) = ii(1)
            jj(7) = jj(3)
            kk(7) = kk(5)

            ! (i+1, j+1, k+1)
            ii(8) = ii(2)
            jj(8) = jj(3)
            kk(8) = kk(5)

        end subroutine trilinear

        pure subroutine trilinear_weights_add(pos, i, j, k, ww)
            double precision, intent(in)    :: pos(3)
            double precision, intent(inout) :: ww(ngp)
            integer,          intent(in)    :: i, j, k
            double precision                :: xyz(3)
            xyz = (pos - lower) * dxi
            call get_weights(xyz, i, j, k, ww)

        end subroutine trilinear_weights_add


        pure subroutine get_weights(xyz, i, j, k, ww)
            double precision, intent(in)    :: xyz(3)
            double precision, intent(inout) :: ww(ngp)
            integer,          intent(in)    :: i, j, k
            double precision                :: px, py, pz, pxc, pyc, pzc
            double precision                :: w00, w10, w01, w11

            ! (i, j, k)
            px = xyz(1) - dble(i)
            pxc = one - px

            py = xyz(2) - dble(j)
            pyc = one - py

            pz = xyz(3) - dble(k)
            pzc = one - pz

            w00 = pyc * pxc
            w10 = pyc * px
            w01 = py * pxc
            w11 = py * px

            ww(1) = ww(1) + pzc * w00  !w000
            ww(2) = ww(2) + pzc * w10  !w100
            ww(3) = ww(3) + pzc * w01  !w010
            ww(4) = ww(4) + pzc * w11  !w110
            ww(5) = ww(5) + pz  * w00  !w001
            ww(6) = ww(6) + pz  * w10  !w101
            ww(7) = ww(7) + pz  * w01  !w011
            ww(8) = ww(8) + pz  * w11  !w111

        end subroutine get_weights

end module parcel_interpl
