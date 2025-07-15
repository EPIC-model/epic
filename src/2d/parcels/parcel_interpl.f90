! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : zero, one, two
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin
    use options, only : parcel
    use dynamic_parcels, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_types, only : idealised_parcel_type, realistic_parcel_type
    use parcel_ellipse
    use fields
    use physics, only : glat, lambda_c, q_0
    use omp_lib
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js

    ! interpolation weights
    double precision :: weights(0:1,0:1)

    integer :: par2grid_timer, &
#ifndef NDBEBUG
               sym_vol2grid_timer, &
#endif
               grid2par_timer

    private :: is, js, weights

    contains

        ! Interpolate the parcel volume to the grid
        subroutine vol2grid
            double precision  :: points(2, 2)
            integer           :: n, p, l
            double precision  :: pvol

            volg = zero

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, is, js, weights) &
            !$omp& reduction(+: volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(:, n), &
                                            pvol, parcels%B(:, n))


                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    volg(js:js+1, is:is+1) = volg(js:js+1, is:is+1) &
                                           + f12 * weights * pvol
                enddo
            enddo
            !$omp end do
            !$omp end parallel
            ! apply periodicity
            volg(:, 0)    = volg(:, 0) + volg(:, nx)
            volg(:, nx-1) = volg(:, nx-1) + volg(:, -1)
            volg(:, -1)   = volg(:, nx-1)
            volg(:, nx)   = volg(:, 0)

            ! apply free slip boundary condition
            volg(0,  :) = two * volg(0,  :)
            volg(nz, :) = two * volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :) = volg(1,    :) + volg(-1,   :)
            volg(nz-1, :) = volg(nz-1, :) + volg(nz+1, :)

        end subroutine vol2grid

#ifndef NDEBUG
        ! Interpolate the parcel volume to the grid to check symmetry
        subroutine vol2grid_symmetry_error
            double precision :: points(2, 2), V, B(2), pos(2)
            integer          :: n, p, l, m
            double precision :: pvol

            call start_timer(sym_vol2grid_timer)

            sym_volg = zero

            do m = -1, 1, 2
                !$omp parallel default(shared)
                !$omp do private(n, p, l, points, pos, pvol, V, B, is, js, weights) &
                !$omp& reduction(+: sym_volg)
                do n = 1, n_parcels

                    pos = parcels%position(:, n)
                    pvol = parcels%volume(n)
                    pos(1) = dble(m) * pos(1)
                    V = dble(m) * pvol
                    B = parcels%B(:, n)

                    B(2) = dble(m) * B(2)

                    points = get_ellipse_points(pos, V, B)

                    ! we have 2 points per ellipse
                    do p = 1, 2

                        ! ensure point is within the domain
                        call apply_periodic_bc(points(:, p))

                        ! get interpolation weights and mesh indices
                        call bilinear(points(:, p), is, js, weights)

                        sym_volg(js:js+1, is:is+1) = sym_volg(js:js+1, is:is+1) &
                                               + f12 * weights * pvol
                    enddo
                enddo
                !$omp end do
                !$omp end parallel
            enddo

            ! apply periodicity
            sym_volg(:, 0)    = sym_volg(:, 0) + sym_volg(:, nx)
            sym_volg(:, nx-1) = sym_volg(:, nx-1) + sym_volg(:, -1)
            sym_volg(:, -1)   = sym_volg(:, nx-1)
            sym_volg(:, nx)   = sym_volg(:, 0)

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
        subroutine par2grid_idealised(parcels)
            class(idealised_parcel_type), intent(in) :: parcels
            double precision :: points(2, 2)
            integer          :: n, p, i, j
            double precision :: pvol, weight(0:1, 0:1), btot

            call start_timer(par2grid_timer)

            vortg = zero
            volg = zero
            nparg = zero
            nsparg = zero
            if(parcels%is_moist) then
                dbuoyg = zero
                humg = zero
            endif
            tbuoyg = zero
            !$omp parallel default(shared)
            !$omp do private(n, p, i, j, points, pvol, weight, btot, is, js, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, dbuoyg, humg, tbuoyg, volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                call parcels%get_buoyancy(n, btot)

                points = get_ellipse_points(parcels%position(:, n), &
                                            pvol, parcels%B(:, n))

                call get_index(parcels%position(:, n), i, j)
                i = mod(i + nx, nx)
                nparg(j, i) = nparg(j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(j, i) = nsparg(j, i) + 1
                endif

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    weight = f12 * weights * pvol

                    vortg(js:js+1, is:is+1) = vortg(js:js+1, is:is+1) &
                                        + weight * parcels%vorticity(1, n)

                    if(parcels%is_moist) then
                        dbuoyg(js:js+1, is:is+1) = dbuoyg(js:js+1, is:is+1) &
                                             + weight * parcels%buoyancy(n)
                        humg(js:js+1, is:is+1) = humg(js:js+1, is:is+1) &
                                           + weight * parcels%humidity(n)
                    endif

                    tbuoyg(js:js+1, is:is+1) = tbuoyg(js:js+1, is:is+1) &
                                         + weight * btot
                    volg(js:js+1, is:is+1) = volg(js:js+1, is:is+1) &
                                       + weight
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ! apply periodicity
            volg(:, 0)    = volg(:, 0) + volg(:, nx)
            volg(:, nx-1) = volg(:, nx-1) + volg(:, -1)
            volg(:, -1)   = volg(:, nx-1)
            volg(:, nx)   = volg(:, 0)

            nparg(:, 0)    = nparg(:, 0) + nparg(:, nx)
            nparg(:, nx-1) = nparg(:, nx-1) + nparg(:, -1)

            nsparg(:, 0)    = nsparg(:, 0) + nsparg(:, nx)
            nsparg(:, nx-1) = nsparg(:, nx-1) + nsparg(:, -1)

            vortg(:, 0)    = vortg(:, 0) + vortg(:, nx)
            vortg(:, nx-1) = vortg(:, nx-1) + vortg(:, -1)
            vortg(:, -1)   = vortg(:, nx-1)
            vortg(:, nx)   = vortg(:, 0)

            tbuoyg(:, 0)    = tbuoyg(:, 0) + tbuoyg(:, nx)
            tbuoyg(:, nx-1) = tbuoyg(:, nx-1) + tbuoyg(:, -1)
            tbuoyg(:, -1)   = tbuoyg(:, nx-1)
            tbuoyg(:, nx)   = tbuoyg(:, 0)

            if(parcels%is_moist) then
                dbuoyg(:, 0)    = dbuoyg(:, 0) + dbuoyg(:, nx)
                dbuoyg(:, nx-1) = dbuoyg(:, nx-1) + dbuoyg(:, -1)
                dbuoyg(:, -1)   = dbuoyg(:, nx-1)
                dbuoyg(:, nx)   = dbuoyg(:, 0)

                humg(:, 0)    = humg(:, 0) + humg(:, nx)
                humg(:, nx-1) = humg(:, nx-1) + humg(:, -1)
                humg(:, -1)   = humg(:, nx-1)
                humg(:, nx)   = humg(:, 0)
            endif

            ! apply free slip boundary condition
            volg(0,  :) = two * volg(0,  :)
            volg(nz, :) = two * volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :) = volg(1,    :) + volg(-1,   :)
            volg(nz-1, :) = volg(nz-1, :) + volg(nz+1, :)

            vortg(0,  :) = two * vortg(0,  :)
            vortg(nz, :) = two * vortg(nz, :)
            vortg(1,    :) = vortg(1,    :) + vortg(-1,   :)
            vortg(nz-1, :) = vortg(nz-1, :) + vortg(nz+1, :)

            if(parcels%is_moist) then
                dbuoyg(0,  :) = two * dbuoyg(0,  :)
                dbuoyg(nz, :) = two * dbuoyg(nz, :)
                dbuoyg(1,    :) = dbuoyg(1,    :) + dbuoyg(-1,   :)
                dbuoyg(nz-1, :) = dbuoyg(nz-1, :) + dbuoyg(nz+1, :)
                humg(0,  :) = two * humg(0,  :)
                humg(nz, :) = two * humg(nz, :)
                humg(1,    :) = humg(1,    :) + humg(-1,   :)
                humg(nz-1, :) = humg(nz-1, :) + humg(nz+1, :)
            endif

            tbuoyg(0,  :) = two * tbuoyg(0,  :)
            tbuoyg(nz, :) = two * tbuoyg(nz, :)
            tbuoyg(1,    :) = tbuoyg(1,    :) + tbuoyg(-1,   :)
            tbuoyg(nz-1, :) = tbuoyg(nz-1, :) + tbuoyg(nz+1, :)
            ! exclude halo cells to avoid division by zero
            vortg(0:nz, :) = vortg(0:nz, :) / volg(0:nz, :)

            ! extrapolate to halo grid points (since halo grid points
            ! are used to get u_z = w_x - zeta)
            vortg(-1,   :) = two * vortg(0,  :) - vortg(1,    :)
            vortg(nz+1, :) = two * vortg(nz, :) - vortg(nz-1, :)

            if(parcels%is_moist) then
                dbuoyg(0:nz, :) = dbuoyg(0:nz, :) / volg(0:nz, :)
                humg(0:nz, :) = humg(0:nz, :) / volg(0:nz, :)
            endif
            tbuoyg(0:nz, :) = tbuoyg(0:nz, :) / volg(0:nz, :)

            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            tbuoyg(-1,   :) = two * tbuoyg(0,  :) - tbuoyg(1, :)
            tbuoyg(nz+1, :) = two * tbuoyg(nz, :) - tbuoyg(nz-1, :)

            ! sum halo contribution into internal cells
            ! (be aware that halo cell contribution at upper boundary
            ! are added to cell nz)
            nparg(0,    :) = nparg(0,    :) + nparg(-1, :)
            nparg(nz-1, :) = nparg(nz-1, :) + nparg(nz, :)

            nsparg(0,    :) = nsparg(0,    :) + nsparg(-1, :)
            nsparg(nz-1, :) = nsparg(nz-1, :) + nsparg(nz, :)

            ! sanity check
            if (sum(nparg(0:nz-1, :)) /= n_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

            call stop_timer(par2grid_timer)

        end subroutine par2grid_idealised


        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - vorticity
        !   - buoyancy
        !   - volume
        ! It also updates the scalar fields:
        !   - nparg, that is the number of parcels per grid cell
        !   - nsparg, that is the number of small parcels per grid cell
        subroutine par2grid_realistic(parcels)
            class(realistic_parcel_type), intent(inout) :: parcels
            double precision :: points(2, 2)
            integer          :: n, p, i, j
            double precision :: pvol, weight(0:1, 0:1), btot

            call parcels%saturation_adjustment

            call start_timer(par2grid_timer)
            vortg = zero
            volg = zero
            nparg = zero
            nsparg = zero
            if(parcels%is_moist) then
                qvg = zero
                qlg = zero
            endif
            if(parcels%has_droplets) then
                Nlg = zero
            endif
            thetag = zero
            tbuoyg = zero
            !$omp parallel default(shared)
            !$omp do private(n, p, i, j, points, pvol, weight, btot, is, js, weights) &
            !$omp& reduction(+:nparg, nsparg, vortg, qvg, qlg, tbuoyg, thetag, Nlg, volg)
            do n = 1, n_parcels
                pvol = parcels%volume(n)

                call parcels%get_buoyancy(n, btot)

                points = get_ellipse_points(parcels%position(:, n), &
                                            pvol, parcels%B(:, n))

                call get_index(parcels%position(:, n), i, j)
                i = mod(i + nx, nx)
                nparg(j, i) = nparg(j, i) + 1
                if (parcels%volume(n) <= vmin) then
                    nsparg(j, i) = nsparg(j, i) + 1
                endif

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    weight = f12 * weights * pvol

                    vortg(js:js+1, is:is+1) = vortg(js:js+1, is:is+1) &
                                        + weight * parcels%vorticity(1, n)

                    if(parcels%is_moist) then
                        qvg(js:js+1, is:is+1) = qvg(js:js+1, is:is+1) &
                                             + weight * parcels%qv(n)
                        qlg(js:js+1, is:is+1) = qlg(js:js+1, is:is+1) &
                                           + weight * parcels%ql(n)
                    endif
                    if(parcels%has_droplets) then
                        Nlg(js:js+1, is:is+1) = Nlg(js:js+1, is:is+1) &
                                             + weight * parcels%Nl(n)
                    endif
                    tbuoyg(js:js+1, is:is+1) = tbuoyg(js:js+1, is:is+1) &
                                         + weight * btot
                    thetag(js:js+1, is:is+1) = thetag(js:js+1, is:is+1) &
                                         + weight * parcels%theta(n)
                    volg(js:js+1, is:is+1) = volg(js:js+1, is:is+1) &
                                       + weight
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ! apply periodicity
            volg(:, 0)    = volg(:, 0) + volg(:, nx)
            volg(:, nx-1) = volg(:, nx-1) + volg(:, -1)
            volg(:, -1)   = volg(:, nx-1)
            volg(:, nx)   = volg(:, 0)

            nparg(:, 0)    = nparg(:, 0) + nparg(:, nx)
            nparg(:, nx-1) = nparg(:, nx-1) + nparg(:, -1)

            nsparg(:, 0)    = nsparg(:, 0) + nsparg(:, nx)
            nsparg(:, nx-1) = nsparg(:, nx-1) + nsparg(:, -1)

            vortg(:, 0)    = vortg(:, 0) + vortg(:, nx)
            vortg(:, nx-1) = vortg(:, nx-1) + vortg(:, -1)
            vortg(:, -1)   = vortg(:, nx-1)
            vortg(:, nx)   = vortg(:, 0)

            tbuoyg(:, 0)    = tbuoyg(:, 0) + tbuoyg(:, nx)
            tbuoyg(:, nx-1) = tbuoyg(:, nx-1) + tbuoyg(:, -1)
            tbuoyg(:, -1)   = tbuoyg(:, nx-1)
            tbuoyg(:, nx)   = tbuoyg(:, 0)

            thetag(:, 0)    = thetag(:, 0) + thetag(:, nx)
            thetag(:, nx-1) = thetag(:, nx-1) + thetag(:, -1)
            thetag(:, -1)   = thetag(:, nx-1)
            thetag(:, nx)   = thetag(:, 0)

            if(parcels%is_moist) then
                qvg(:, 0)    = qvg(:, 0) + qvg(:, nx)
                qvg(:, nx-1) = qvg(:, nx-1) + qvg(:, -1)
                qvg(:, -1)   = qvg(:, nx-1)
                qvg(:, nx)   = qvg(:, 0)

                qlg(:, 0)    = qlg(:, 0) + qlg(:, nx)
                qlg(:, nx-1) = qlg(:, nx-1) + qlg(:, -1)
                qlg(:, -1)   = qlg(:, nx-1)
                qlg(:, nx)   = qlg(:, 0)
            endif

            if(parcels%has_droplets) then
                Nlg(:, 0)    = Nlg(:, 0) + Nlg(:, nx)
                Nlg(:, nx-1) = Nlg(:, nx-1) + Nlg(:, -1)
                Nlg(:, -1)   = Nlg(:, nx-1)
                Nlg(:, nx)   = Nlg(:, 0)
            endif

            ! apply free slip boundary condition
            volg(0,  :) = two * volg(0,  :)
            volg(nz, :) = two * volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :) = volg(1,    :) + volg(-1,   :)
            volg(nz-1, :) = volg(nz-1, :) + volg(nz+1, :)

            vortg(0,  :) = two * vortg(0,  :)
            vortg(nz, :) = two * vortg(nz, :)
            vortg(1,    :) = vortg(1,    :) + vortg(-1,   :)
            vortg(nz-1, :) = vortg(nz-1, :) + vortg(nz+1, :)

            if(parcels%is_moist) then
                qvg(0,  :) = two * qvg(0,  :)
                qvg(nz, :) = two * qvg(nz, :)
                qvg(1,    :) = qvg(1,    :) + qvg(-1,   :)
                qvg(nz-1, :) = qvg(nz-1, :) + qvg(nz+1, :)
                qlg(0,  :) = two * qlg(0,  :)
                qlg(nz, :) = two * qlg(nz, :)
                qlg(1,    :) = qlg(1,    :) + qlg(-1,   :)
                qlg(nz-1, :) = qlg(nz-1, :) + qlg(nz+1, :)
            endif

            if(parcels%has_droplets) then
                Nlg(0,  :) = two * Nlg(0,  :)
                Nlg(nz, :) = two * Nlg(nz, :)
                Nlg(1,    :) = Nlg(1,    :) + Nlg(-1,   :)
                Nlg(nz-1, :) = Nlg(nz-1, :) + Nlg(nz+1, :)
            endif

            tbuoyg(0,  :) = two * tbuoyg(0,  :)
            tbuoyg(nz, :) = two * tbuoyg(nz, :)
            tbuoyg(1,    :) = tbuoyg(1,    :) + tbuoyg(-1,   :)
            tbuoyg(nz-1, :) = tbuoyg(nz-1, :) + tbuoyg(nz+1, :)

            thetag(0,  :) = two * thetag(0,  :)
            thetag(nz, :) = two * thetag(nz, :)
            thetag(1,    :) = thetag(1,    :) + thetag(-1,   :)
            thetag(nz-1, :) = thetag(nz-1, :) + thetag(nz+1, :)

            ! exclude halo cells to avoid division by zero
            vortg(0:nz, :) = vortg(0:nz, :) / volg(0:nz, :)

            ! extrapolate to halo grid points (since halo grid points
            ! are used to get u_z = w_x - zeta)
            vortg(-1,   :) = two * vortg(0,  :) - vortg(1,    :)
            vortg(nz+1, :) = two * vortg(nz, :) - vortg(nz-1, :)

            if(parcels%is_moist) then
                qvg(0:nz, :) = qvg(0:nz, :) / volg(0:nz, :)
                qlg(0:nz, :) = qlg(0:nz, :) / volg(0:nz, :)
            endif

            if(parcels%has_droplets) then
                Nlg(0:nz, :) = Nlg(0:nz, :) / volg(0:nz, :)
            endif
            tbuoyg(0:nz, :) = tbuoyg(0:nz, :) / volg(0:nz, :)
            thetag(0:nz, :) = thetag(0:nz, :) / volg(0:nz, :)

            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            tbuoyg(-1,   :) = two * tbuoyg(0,  :) - tbuoyg(1, :)
            tbuoyg(nz+1, :) = two * tbuoyg(nz, :) - tbuoyg(nz-1, :)
            thetag(-1,   :) = two * thetag(0,  :) - thetag(1, :)
            thetag(nz+1, :) = two * thetag(nz, :) - thetag(nz-1, :)

            ! sum halo contribution into internal cells
            ! (be aware that halo cell contribution at upper boundary
            ! are added to cell nz)
            nparg(0,    :) = nparg(0,    :) + nparg(-1, :)
            nparg(nz-1, :) = nparg(nz-1, :) + nparg(nz, :)

            nsparg(0,    :) = nsparg(0,    :) + nsparg(-1, :)
            nsparg(nz-1, :) = nsparg(nz-1, :) + nsparg(nz, :)

            ! sanity check
            if (sum(nparg(0:nz-1, :)) /= n_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

            call stop_timer(par2grid_timer)

        end subroutine par2grid_realistic


        ! Interpolate the gridded quantities to the parcels
        ! @param[inout] vel is the parcel velocity
        ! @param[inout] vor is the parcel vorticity
        ! @param[inout] vgrad is the parcel strain
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine grid2par(vel, vor, vgrad, add)
            double precision,     intent(inout) :: vel(:, :), vor(:, :), vgrad(:, :)
            logical, optional, intent(in)       :: add
            double precision                    :: points(2, 2), weight(0:1, 0:1)
            integer                             :: n, p, l

            call start_timer(grid2par_timer)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_parcels
                        vel(:, n) = zero
                        vor(1, n)    = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    vel(:, n) = zero
                    vor(1, n)    = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, weight, is, js, weights)
            do n = 1, n_parcels

                vgrad(:, n) = zero

                points = get_ellipse_points(parcels%position(:, n), &
                                            parcels%volume(n),      &
                                            parcels%B(:, n))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(:, p))

                    ! get interpolation weights and mesh indices
                    call bilinear(points(:, p), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    weight = f12 * weights

                    ! the weight is halved due to 2 points per ellipse
                    do l = 1,2
                        vel(l, n) = vel(l, n) &
                                  + sum(weight * velog(js:js+1, is:is+1, l))
                    end do
                    do l=1,4
                        vgrad(l, n) = vgrad(l, n) &
                                    + sum(weight * velgradg(js:js+1, is:is+1, l))
                    end do
                    vor(1, n) = vor(1, n) + sum(weight * vtend(js:js+1, is:is+1))
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Bi-linear interpolation
        ! @param[in] pos position vector
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj vertical grid points for interpolation
        ! @param[out] ww interpolation weights
        pure subroutine bilinear(pos, ii, jj, ww)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ii, jj
            double precision, intent(out) :: ww(0:1, 0:1)
            double precision              :: xz(2)
            double precision              :: px, pz, pxc, pzc


            ! (i, j)
            xz = (pos - lower(1:2)) * dxi(1:2)
            ii = floor(xz(1))
            jj = floor(xz(2))

            px = xz(1) - dble(ii)
            pxc = one - px

            pz = xz(2) - dble(jj)
            pzc = one - pz

            ! Note order of indices is j,i
            ww(0, 0) = pzc * pxc
            ww(0, 1) = pzc * px
            ww(1, 0) = pz  * pxc
            ww(1, 1) = pz  * px

        end subroutine bilinear

        ! Bi-linear interpolation
        ! @param[in] pos position of the parcel
        ! @param[out] ii horizontal grid points for interoplation
        ! @param[out] jj vertical grid points for interpolation
        ! @param[out] ww interpolation weights
        subroutine bilinear_old(pos, ii, jj, ww)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ii(4), jj(4)
            double precision, intent(out) :: ww(4)
            double precision              :: xy(2)

            ! (i, j)
            call get_index(pos, ii(1), jj(1))
            call get_position(ii(1), jj(1), xy)
            ww(1) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j)
            ii(2) = ii(1) + 1
            jj(2) = jj(1)
            call get_position(ii(2), jj(2), xy)
            ww(2) = product(one - abs(pos - xy) * dxi)

            ! (i, j+1)
            ii(3) = ii(1)
            jj(3) = jj(1) + 1
            call get_position(ii(3), jj(3), xy)
            ww(3) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j+1)
            ii(4) = ii(2)
            jj(4) = jj(3)
            call get_position(ii(4), jj(4), xy)
            ww(4) = product(one - abs(pos - xy) * dxi)

            ! account for x periodicity
            call periodic_index_shift(ii)

        end subroutine bilinear_old

        subroutine par2grid
            select type (parcels)
            type is (idealised_parcel_type)
                call par2grid_idealised(parcels)
            type is (realistic_parcel_type)
                call par2grid_realistic(parcels)
            end select
        end subroutine par2grid

end module parcel_interpl
