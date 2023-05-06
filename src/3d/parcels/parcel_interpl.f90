! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : zero, one, two, f14
    use mpi_timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin, l_bndry_zeta_zero
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipsoid
    use fields
    use field_mpi, only : field_halo_swap
    use physics, only : glat, lambda_c, q_0
    use omp_lib
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = 8

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(ngp)

    integer :: par2grid_timer, &
               grid2par_timer

    private :: is, js, ks, weights

    contains

        ! Interpolate the parcel volume to the grid
        ! @pre The parcel must be assigned to the correct MPI process.
        subroutine vol2grid(l_reuse)
            logical, optional, intent(in) :: l_reuse
            double precision              :: points(3, 4)
            integer                       :: n, p
            double precision              :: pvol, f14pvol

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, points, pvol, is, js, ks, weights, f14pvol) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n),  &
                                              n, l_reuse)

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    call trilinear(points(:, p), is, js, ks, weights)
                    
                    f14pvol = f14 * pvol

                    volg(ks, js, is) = volg(ks, js, is) &
                                                  + f14pvol * weights(1)
                    volg(ks, js, is+1) = volg(ks, js, is+1) &
                                                  + f14pvol * weights(2)
                    volg(ks, js+1, is) = volg(ks, js+1, is) &
                                                  + f14pvol * weights(3)
                    volg(ks, js+1, is+1) = volg(ks, js+1, is+1) &
                                                  + f14pvol * weights(4)
                    volg(ks+1, js, is) = volg(ks+1, js, is) &
                                                  + f14pvol * weights(5)
                    volg(ks+1, js, is+1) = volg(ks+1, js, is+1) &
                                                  + f14pvol * weights(6)
                    volg(ks+1, js+1, is) = volg(ks+1, js+1, is) &
                                                  + f14pvol * weights(7)
                    volg(ks+1, js+1, is+1) = volg(ks+1, js+1, is+1) &
                                                  + f14pvol * weights(8)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call field_halo_swap(volg)

            ! apply free slip boundary condition
            !$omp parallel workshare
            volg(0,  :, :) = two * volg(0,  :, :)
            volg(nz, :, :) = two * volg(nz, :, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :, :) = volg(1,    :, :) + volg(-1,   :, :)
            volg(nz-1, :, :) = volg(nz-1, :, :) + volg(nz+1, :, :)
            !$omp end parallel workshare

        end subroutine vol2grid


        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        !   - volume
        ! It also updates the scalar fields:
        !   - nparg, that is the number of parcels per grid cell
        !   - nsparg, that is the number of small parcels per grid cell
        ! @pre The parcel must be assigned to the correct MPI process.
        subroutine par2grid(l_reuse)
            logical, optional :: l_reuse
            double precision  :: points(3, 4)
            integer           :: n, p, i, j, k
            double precision  :: pvol, btot, f14pvol
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
            !$omp do private(n, p, i, j, k, points, pvol, btot, q_c, f14pvol) &
            !$omp& private( is, js, ks, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, dbuoyg, humg, tbuoyg, volg)
#else
            !$omp do private(n, p, i, j, k, points, pvol, btot, f14pvol) &
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
                nparg(k, j, i) = nparg(k, j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(k, j, i) = nsparg(k, j, i) + 1
                endif

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is a quarter due to 4 points per ellipsoid

                    f14pvol = f14 * pvol

                    vortg(ks, js, is, :) = vortg(ks, js, is, :) &
                                                  + f14pvol * weights(1) * parcels%vorticity(:, n)
                    vortg(ks, js, is+1, :) = vortg(ks, js, is+1, :) &
                                                  + f14pvol * weights(2) * parcels%vorticity(:, n)
                    vortg(ks, js+1, is, :) = vortg(ks, js+1, is, :) &
                                                  + f14pvol * weights(3) * parcels%vorticity(:, n)
                    vortg(ks, js+1, is+1, :) = vortg(ks, js+1, is+1, :) &
                                                  + f14pvol * weights(4) * parcels%vorticity(:, n)
                    vortg(ks+1, js, is, :) = vortg(ks+1, js, is, :) &
                                                  + f14pvol * weights(5) * parcels%vorticity(:, n)
                    vortg(ks+1, js, is+1, :) = vortg(ks+1, js, is+1, :) &
                                                  + f14pvol * weights(6) * parcels%vorticity(:, n)
                    vortg(ks+1, js+1, is, :) = vortg(ks+1, js+1, is, :) &
                                                  + f14pvol * weights(7) * parcels%vorticity(:, n)
                    vortg(ks+1, js+1, is+1, :) = vortg(ks+1, js+1, is+1, :) &
                                                  + f14pvol * weights(8) * parcels%vorticity(:, n)

#ifndef ENABLE_DRY_MODE
                    dbuoyg(ks, js, is) = dbuoyg(ks, js, is) &
                                                  + f14pvol * weights(1) * parcels%buoyancy(n)
                    dbuoyg(ks, js, is+1) = dbuoyg(ks, js, is+1) &
                                                  + f14pvol * weights(2) * parcels%buoyancy(n)
                    dbuoyg(ks, js+1, is) = dbuoyg(ks, js+1, is) &
                                                  + f14pvol * weights(3) * parcels%buoyancy(n)
                    dbuoyg(ks, js+1, is+1) = dbuoyg(ks, js+1, is+1) &
                                                  + f14pvol * weights(4) * parcels%buoyancy(n)
                    dbuoyg(ks+1, js, is) = dbuoyg(ks+1, js, is) &
                                                  + f14pvol * weights(5) * parcels%buoyancy(n)
                    dbuoyg(ks+1, js, is+1) = dbuoyg(ks+1, js, is+1) &
                                                  + f14pvol * weights(6) * parcels%buoyancy(n)
                    dbuoyg(ks+1, js+1, is) = dbuoyg(ks+1, js+1, is) &
                                                  + f14pvol * weights(7) * parcels%buoyancy(n)
                    dbuoyg(ks+1, js+1, is+1) = dbuoyg(ks+1, js+1, is+1) &
                                                  + f14pvol * weights(8) * parcels%buoyancy(n)

                    humg(ks, js, is) = humg(ks, js, is) &
                                                  + f14pvol * weights(1) * parcels%humidity(n)
                    humg(ks, js, is+1) = humg(ks, js, is+1) &
                                                  + f14pvol * weights(2) * parcels%humidity(n)
                    humg(ks, js+1, is) = humg(ks, js+1, is) &
                                                  + f14pvol * weights(3) * parcels%humidity(n)
                    humg(ks, js+1, is+1) = humg(ks, js+1, is+1) &
                                                  + f14pvol * weights(4) * parcels%humidity(n)
                    humg(ks+1, js, is) = humg(ks+1, js, is) &
                                                  + f14pvol * weights(5) * parcels%humidity(n)
                    humg(ks+1, js, is+1) = humg(ks+1, js, is+1) &
                                                  + f14pvol * weights(6) * parcels%humidity(n)
                    humg(ks+1, js+1, is) = humg(ks+1, js+1, is) &
                                                  + f14pvol * weights(7) * parcels%humidity(n)
                    humg(ks+1, js+1, is+1) = humg(ks+1, js+1, is+1) &
                                                  + f14pvol * weights(8) * parcels%humidity(n)
#endif
                    tbuoyg(ks, js, is) = tbuoyg(ks, js, is) &
                                                  + f14pvol * weights(1) * btot
                    tbuoyg(ks, js, is+1) = tbuoyg(ks, js, is+1) &
                                                  + f14pvol * weights(2) * btot
                    tbuoyg(ks, js+1, is) = tbuoyg(ks, js+1, is) &
                                                  + f14pvol * weights(3) * btot
                    tbuoyg(ks, js+1, is+1) = tbuoyg(ks, js+1, is+1) &
                                                  + f14pvol * weights(4) * btot
                    tbuoyg(ks+1, js, is) = tbuoyg(ks+1, js, is) &
                                                  + f14pvol * weights(5) * btot
                    tbuoyg(ks+1, js, is+1) = tbuoyg(ks+1, js, is+1) &
                                                  + f14pvol * weights(6) * btot
                    tbuoyg(ks+1, js+1, is) = tbuoyg(ks+1, js+1, is) &
                                                  + f14pvol * weights(7) * btot
                    tbuoyg(ks+1, js+1, is+1) = tbuoyg(ks+1, js+1, is+1) &
                                                  + f14pvol * weights(8) * btot

                    volg(ks, js, is) = volg(ks, js, is) &
                                                  + f14pvol * weights(1)
                    volg(ks, js, is+1) = volg(ks, js, is+1) &
                                                  + f14pvol * weights(2)
                    volg(ks, js+1, is) = volg(ks, js+1, is) &
                                                  + f14pvol * weights(3)
                    volg(ks, js+1, is+1) = volg(ks, js+1, is+1) &
                                                  + f14pvol * weights(4)
                    volg(ks+1, js, is) = volg(ks+1, js, is) &
                                                  + f14pvol * weights(5)
                    volg(ks+1, js, is+1) = volg(ks+1, js, is+1) &
                                                  + f14pvol * weights(6)
                    volg(ks+1, js+1, is) = volg(ks+1, js+1, is) &
                                                  + f14pvol * weights(7)
                    volg(ks+1, js+1, is+1) = volg(ks+1, js+1, is+1) &
                                                  + f14pvol * weights(8)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call field_halo_swap(volg)
            call field_halo_swap(vortg(:, :, :, I_X))
            call field_halo_swap(vortg(:, :, :, I_Y))
            call field_halo_swap(vortg(:, :, :, I_Z))
            call field_halo_swap(tbuoyg)

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
            vortg(1,    :, :, :) = vortg(1,    :, :, :) + vortg(-1,   :, :, :)
            vortg(nz-1, :, :, :) = vortg(nz-1, :, :, :) + vortg(nz+1, :, :, :)

            tbuoyg(0,  :, :) = two * tbuoyg(0,  :, :)
            tbuoyg(nz, :, :) = two * tbuoyg(nz, :, :)
            tbuoyg(1,    :, :) = tbuoyg(1,    :, :) + tbuoyg(-1,   :, :)
            tbuoyg(nz-1, :, :) = tbuoyg(nz-1, :, :) + tbuoyg(nz+1, :, :)
            !$omp end parallel workshare

#ifndef ENABLE_DRY_MODE
            call field_halo_swap(dbuoyg)
            call field_halo_swap(humg)
            !$omp parallel workshare
            dbuoyg(0,  :, :) = two * dbuoyg(0,  :, :)
            dbuoyg(nz, :, :) = two * dbuoyg(nz, :, :)
            dbuoyg(1,    :, :) = dbuoyg(1,    :, :) + dbuoyg(-1,   :, :)
            dbuoyg(nz-1, :, :) = dbuoyg(nz-1, :, :) + dbuoyg(nz+1, :, :)
            humg(0,  :, :) = two * humg(0,  :, :)
            humg(nz, :, :) = two * humg(nz, :, :)
            humg(1,    :, :) = humg(1,    :, :) + humg(-1,   :, :)
            humg(nz-1, :, :) = humg(nz-1, :, :) + humg(nz+1, :, :)
            !$omp end parallel workshare
#endif

            ! exclude halo cells to avoid division by zero
            do p = 1, 3
                vortg(0:nz, :, :, p) = vortg(0:nz, :, :, p) / volg(0:nz, :, :)
            enddo

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
                call mpi_exit_on_error("par2grid: Wrong total number of parcels!")
            endif

            call stop_timer(par2grid_timer)

        end subroutine par2grid


        ! Interpolate the gridded quantities to the parcels
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        ! @pre The parcel must be assigned to the correct MPI process and the halo of fields must be
        !      filled correctly.
        subroutine grid2par(add)
            logical, optional, intent(in) :: add
            double precision              :: points(3, 4)
            integer                       :: n, p

            call start_timer(grid2par_timer)

            ! clear old data efficiently
            if (present(add)) then
               if (add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_parcels
                        parcels%delta_pos(:, n) = zero
                        parcels%delta_vor(:, n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    parcels%delta_pos(:, n) = zero
                    parcels%delta_vor(:, n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, points, is, js, ks, weights)
            do n = 1, n_parcels

                parcels%strain(:, n) = zero

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              parcels%volume(n), parcels%B(:, n), n)

                do p = 1, 4
                    ! get interpolation weights and mesh indices
                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    parcels%delta_pos(:, n) = parcels%delta_pos(:, n) &
                                                + f14 * (weights(1) * velog(ks, js, is, :) &
                                                + weights(2) * velog(ks  , js   , is+1, :) &
                                                + weights(3) * velog(ks  , js+1 , is,   :) &
                                                + weights(4) * velog(ks  , js+1 , is+1, :) &
                                                + weights(5) * velog(ks+1, js   , is,   :) &
                                                + weights(6) * velog(ks+1, js   , is,   :) &
                                                + weights(7) * velog(ks+1, js+1 , is,   :) &
                                                + weights(8) * velog(ks+1, js+1 , is+1, :) )

                    parcels%strain(:, n) = parcels%strain(:, n) &
                                             + f14 * (weights(1) * velgradg(ks, js, is, :) &
                                                + weights(2) * velgradg(ks  , js   , is+1, :) &
                                                + weights(3) * velgradg(ks  , js+1 , is,   :) &
                                                + weights(4) * velgradg(ks  , js+1 , is+1, :) &
                                                + weights(5) * velgradg(ks+1, js   , is,   :) &
                                                + weights(6) * velgradg(ks+1, js   , is,   :) &
                                                + weights(7) * velgradg(ks+1, js+1 , is,   :) &
                                                + weights(8) * velgradg(ks+1, js+1 , is+1, :) )

                    parcels%delta_vor(:, n) = parcels%delta_vor(:, n) &
                                                + f14 * (weights(1) * vtend(ks, js, is, :) &
                                                + weights(2) * vtend(ks  , js   , is+1, :) &
                                                + weights(3) * vtend(ks  , js+1 , is,   :) &
                                                + weights(4) * vtend(ks  , js+1 , is+1, :) &
                                                + weights(5) * vtend(ks+1, js   , is,   :) &
                                                + weights(6) * vtend(ks+1, js   , is,   :) &
                                                + weights(7) * vtend(ks+1, js+1 , is,   :) &
                                                + weights(8) * vtend(ks+1, js+1 , is+1, :) )

                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(grid2par_timer)

        end subroutine grid2par

        ! Tri-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii zonal grid points for interoplation
        ! @param[out] jj meridional grid points for interpolation
        ! @param[out] kk vertical grid points for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine trilinear(pos, ii, jj, kk, ww)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: ii, jj, kk
            double precision, intent(out) :: ww(ngp)
            double precision              :: xyz(3)

            ww = zero

            ! (i, j, k)
            xyz = (pos - lower) * dxi
            ii = floor(xyz(1))
            jj = floor(xyz(2))
            kk = floor(xyz(3))

            call get_weights(xyz, ii, jj, kk, ww)

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
