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
    use field_mpi, only : field_mpi_alloc                   &
                        , field_mpi_dealloc                 &
                        , field_buffer_to_halo              &
                        , field_halo_to_buffer              &
                        , field_buffer_to_interior          &
                        , field_interior_to_buffer          &
                        , interior_to_halo_communication    &
                        , halo_to_interior_communication    &
                        , field_halo_swap
    use physics, only : glat, lambda_c, q_0
    use omp_lib
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks

    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)

    integer :: par2grid_timer, &
               grid2par_timer, &
               halo_swap_timer

    integer, parameter :: IDX_VOL_SWAP   = 1    &
                        , IDX_VOR_X_SWAP = 2    &
                        , IDX_VOR_Y_SWAP = 3    &
                        , IDX_VOR_Z_SWAP = 4    &
                        , IDX_TBUOY_SWAP = 5
#ifndef ENABLE_DRY_MODE
    integer, parameter :: IDX_DBUOY_SWAP = 6    &
                        , IDX_HUM_SWAP   = 7

    integer, parameter :: n_field_swap = 7
#else
    integer, parameter :: n_field_swap = 5
#endif

    public :: par2grid          &
            , vol2grid          &
            , grid2par          &
            , par2grid_timer    &
            , grid2par_timer    &
            , halo_swap_timer    &
            , trilinear         

    contains

        ! Interpolate the parcel volume to the grid
        ! @pre The parcel must be assigned to the correct MPI process.
        subroutine vol2grid(l_reuse)
            logical, optional, intent(in) :: l_reuse
            double precision              :: points(3, 4)
            integer                       :: n, p
            double precision              :: pvol

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, points, pvol, is, js, ks, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              pvol, parcels%B(:, n),  &
                                              n, l_reuse)

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    call trilinear(points(:, p), is, js, ks, weights)

                    volg(ks:ks+1, js:js+1, is:is+1) = volg(ks:ks+1, js:js+1, is:is+1) &
                                              + f14 * weights * pvol
                enddo
            enddo
            !$omp end do
            !$omp end parallel
            
            call start_timer(halo_swap_timer)
            call field_halo_swap(volg)
            call stop_timer(halo_swap_timer)

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
            integer           :: n, p, l, i, j, k
            double precision  :: pvol, weight(0:1,0:1,0:1), btot
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
                nparg(k, j, i) = nparg(k, j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(k, j, i) = nsparg(k, j, i) + 1
                endif

                ! we have 4 points per ellipsoid
                do p = 1, 4

                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is a quarter due to 4 points per ellipsoid
                    weight = f14 * pvol* weights

                    do l = 1, 3
                        vortg(ks:ks+1, js:js+1, is:is+1, l) = vortg(ks:ks+1, js:js+1, is:is+1, l) &
                                                            + weight * parcels%vorticity(l, n)
                    enddo
#ifndef ENABLE_DRY_MODE
                    dbuoyg(ks:ks+1, js:js+1, is:is+1) = dbuoyg(ks:ks+1, js:js+1, is:is+1) &
                                                      + weight * parcels%buoyancy(n)
                    humg(ks:ks+1, js:js+1, is:is+1) = humg(ks:ks+1, js:js+1, is:is+1) &
                                              + weight * parcels%humidity(n)
#endif
                    tbuoyg(ks:ks+1, js:js+1, is:is+1) = tbuoyg(ks:ks+1, js:js+1, is:is+1) &
                                                + weight * btot
                    volg(ks:ks+1, js:js+1, is:is+1) = volg(ks:ks+1, js:js+1, is:is+1) &
                                              + weight
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call start_timer(halo_swap_timer)
            call par2grid_halo_swap
            call stop_timer(halo_swap_timer)

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


        subroutine par2grid_halo_swap
            ! we must first fill the interior grid points
            ! correctly, and then the halo; otherwise
            ! halo grid points do not have correct values at
            ! corners where multiple processes share grid points.

            call field_mpi_alloc(n_field_swap)

            !------------------------------------------------------------------
            ! Accumulate interior:

            call field_halo_to_buffer(volg,                IDX_VOL_SWAP)
            call field_halo_to_buffer(vortg(:, :, :, I_X), IDX_VOR_X_SWAP)
            call field_halo_to_buffer(vortg(:, :, :, I_Y), IDX_VOR_Y_SWAP)
            call field_halo_to_buffer(vortg(:, :, :, I_Z), IDX_VOR_Z_SWAP)
            call field_halo_to_buffer(tbuoyg,              IDX_TBUOY_SWAP)
#ifndef ENABLE_DRY_MODE
            call field_halo_to_buffer(dbuoyg,              IDX_DBUOY_SWAP)
            call field_halo_to_buffer(humg,                IDX_HUM_SWAP)
#endif

            ! send halo data to valid regions of other processes
            call halo_to_interior_communication

            ! accumulate interior; after this operation
            ! all interior grid points have the correct value
            call field_buffer_to_interior(volg,                IDX_VOL_SWAP, .true.)
            call field_buffer_to_interior(vortg(:, :, :, I_X), IDX_VOR_X_SWAP, .true.)
            call field_buffer_to_interior(vortg(:, :, :, I_Y), IDX_VOR_Y_SWAP, .true.)
            call field_buffer_to_interior(vortg(:, :, :, I_Z), IDX_VOR_Z_SWAP, .true.)
            call field_buffer_to_interior(tbuoyg,              IDX_TBUOY_SWAP, .true.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_interior(dbuoyg,              IDX_DBUOY_SWAP, .true.)
            call field_buffer_to_interior(humg,                IDX_HUM_SWAP, .true.)
#endif

            !------------------------------------------------------------------
            ! Fill halo:

            call field_interior_to_buffer(volg,                IDX_VOL_SWAP)
            call field_interior_to_buffer(vortg(:, :, :, I_X), IDX_VOR_X_SWAP)
            call field_interior_to_buffer(vortg(:, :, :, I_Y), IDX_VOR_Y_SWAP)
            call field_interior_to_buffer(vortg(:, :, :, I_Z), IDX_VOR_Z_SWAP)
            call field_interior_to_buffer(tbuoyg,              IDX_TBUOY_SWAP)
#ifndef ENABLE_DRY_MODE
            call field_interior_to_buffer(dbuoyg,              IDX_DBUOY_SWAP)
            call field_interior_to_buffer(humg,                IDX_HUM_SWAP)
#endif

            call interior_to_halo_communication

            call field_buffer_to_halo(volg,                IDX_VOL_SWAP, .false.)
            call field_buffer_to_halo(vortg(:, :, :, I_X), IDX_VOR_X_SWAP, .false.)
            call field_buffer_to_halo(vortg(:, :, :, I_Y), IDX_VOR_Y_SWAP, .false.)
            call field_buffer_to_halo(vortg(:, :, :, I_Z), IDX_VOR_Z_SWAP, .false.)
            call field_buffer_to_halo(tbuoyg,              IDX_TBUOY_SWAP, .false.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_halo(dbuoyg,              IDX_DBUOY_SWAP, .false.)
            call field_buffer_to_halo(humg,                IDX_HUM_SWAP, .false.)
#endif

            call field_mpi_dealloc

        end subroutine par2grid_halo_swap


        ! Interpolate the gridded quantities to the parcels
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        ! @pre The parcel must be assigned to the correct MPI process and the halo of fields must be
        !      filled correctly.
        subroutine grid2par(add)
            logical, optional, intent(in) :: add
            double precision              :: points(3, 4)
            integer                       :: n, l, p

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
            !$omp do private(n, l, p, points, is, js, ks, weights)
            do n = 1, n_parcels

                parcels%strain(:, n) = zero

                points = get_ellipsoid_points(parcels%position(:, n), &
                                              parcels%volume(n), parcels%B(:, n), n)

                do p = 1, 4
                    ! get interpolation weights and mesh indices
                    call trilinear(points(:, p), is, js, ks, weights)

                    ! loop over grid points which are part of the interpolation
                    do l = 1,3
                        parcels%delta_pos(l, n) = parcels%delta_pos(l, n) &
                                                + f14 * sum(weights * velog(ks:ks+1, js:js+1, is:is+1, l))
                    enddo
                    do l = 1,5
                        parcels%strain(l, n) = parcels%strain(l, n) &
                                             + f14 * sum(weights * velgradg(ks:ks+1, js:js+1, is:is+1, l))
                    enddo
                    do l = 1,3
                        parcels%delta_vor(l, n) = parcels%delta_vor(l, n) &
                                                + f14 * sum(weights * vtend(ks:ks+1, js:js+1, is:is+1, l))
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(grid2par_timer)

        end subroutine grid2par

        ! Tri-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii zonal lower grid point for interoplation
        ! @param[out] jj meridional lower grid point for interpolation
        ! @param[out] kk vertical lower grid point for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine trilinear(pos, ii, jj, kk, ww)
            double precision, intent(in)  :: pos(3)
            integer,          intent(out) :: ii, jj, kk
            double precision, intent(out) :: ww(0:1,0:1,0:1)
            double precision              :: xyz(3)

            ww = zero

            ! (i, j, k)
            xyz = (pos - lower) * dxi
            ii = floor(xyz(1))
            jj = floor(xyz(2))
            kk = floor(xyz(3))

            call get_weights(xyz, ii, jj, kk, ww)

        end subroutine trilinear


        pure subroutine get_weights(xyz, i, j, k, ww)
            double precision, intent(in)    :: xyz(3)
            double precision, intent(inout) :: ww(0:1,0:1,0:1)
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

            ! Note order of indices is k,j,i
            ww(0,0,0) = ww(0,0,0) + pzc * w00
            ww(0,0,1) = ww(0,0,1) + pzc * w10
            ww(0,1,0) = ww(0,1,0) + pzc * w01
            ww(0,1,1) = ww(0,1,1) + pzc * w11
            ww(1,0,0) = ww(1,0,0) + pz * w00
            ww(1,0,1) = ww(1,0,1) + pz * w10
            ww(1,1,0) = ww(1,1,0) + pz * w01
            ww(1,1,1) = ww(1,1,1) + pz * w11

        end subroutine get_weights


end module parcel_interpl
