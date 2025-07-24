! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module prec_parcel_interpl
    use constants, only : zero, two
    use timer, only : start_timer, stop_timer
    use parameters, only : nx, nz
    use precipitation_parcels, only : prec_parcels, n_prec_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_types, only : prec_parcel_type
    use parcel_interpl, only :  bilinear
    use fields
    use omp_lib
    implicit none

    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js

    ! interpolation weights
    double precision :: weights(0:1,0:1)

    integer :: prec_par2grid_timer, &
               prec_grid2par_timer

    private :: is, js, weights

    contains

        ! Interpolate parcel quantities to the grid, these consist of the parcel
        !   - qr
        !   - Nr
        !   - precipitation volume (just a diagnostic really)
        ! It also updates the scalar fields:
        !   - prec_nparg, that is the number of parcels per grid cell
        subroutine prec_par2grid(prec_parcels)
            class(prec_parcel_type), intent(in) :: prec_parcels
            double precision :: points(2)
            integer          :: n, p, i, j
            double precision :: pvol, btot

            call start_timer(prec_par2grid_timer)

            prec_volg = zero
            prec_nparg = zero
            prec_tbuoyg = zero
            qrg = zero
            Nrg = zero
            !$omp parallel default(shared)
            !$omp do private(n, p, i, j, points, pvol, btot, is, js, weights) &
            !$omp& reduction(+:prec_nparg, qrg, Nrg, prec_tbuoyg, prec_volg)
            do n = 1, n_prec_parcels
                pvol = prec_parcels%volume(n)

                call prec_parcels%get_buoyancy(n, btot)

                points = prec_parcels%position(:, n)

                call get_index(prec_parcels%position(:, n), i, j)
                i = mod(i + nx, nx)
                prec_nparg(j, i) = prec_nparg(j, i) + 1

                ! ensure point is within the domain
                call apply_periodic_bc(points(:))

                ! get interpolation weights and mesh indices
                call bilinear(points(:), is, js, weights)

                weights = weights*pvol

                prec_tbuoyg(js:js+1, is:is+1) = prec_tbuoyg(js:js+1, is:is+1) &
                                     + weights * btot
                prec_volg(js:js+1, is:is+1) = prec_volg(js:js+1, is:is+1) &
                                   + weights
                qrg(js:js+1, is:is+1) = qrg(js:js+1, is:is+1) &
                                   + weights * prec_parcels%qr(n)
                Nrg(js:js+1, is:is+1) = Nrg(js:js+1, is:is+1) &
                                   + weights * prec_parcels%Nr(n)
            enddo
            !$omp end do
            !$omp end parallel

            ! apply periodicity
            prec_volg(:, 0)    = prec_volg(:, 0) + prec_volg(:, nx)
            prec_volg(:, nx-1) = prec_volg(:, nx-1) + prec_volg(:, -1)
            prec_volg(:, -1)   = prec_volg(:, nx-1)
            prec_volg(:, nx)   = prec_volg(:, 0)

            prec_nparg(:, 0)    = prec_nparg(:, 0) + prec_nparg(:, nx)
            prec_nparg(:, nx-1) = prec_nparg(:, nx-1) + prec_nparg(:, -1)

            prec_tbuoyg(:, 0)    = prec_tbuoyg(:, 0) + prec_tbuoyg(:, nx)
            prec_tbuoyg(:, nx-1) = prec_tbuoyg(:, nx-1) + prec_tbuoyg(:, -1)
            prec_tbuoyg(:, -1)   = prec_tbuoyg(:, nx-1)
            prec_tbuoyg(:, nx)   = prec_tbuoyg(:, 0)

            qrg(:, 0)    = qrg(:, 0) + qrg(:, nx)
            qrg(:, nx-1) = qrg(:, nx-1) + qrg(:, -1)
            qrg(:, -1)   = qrg(:, nx-1)
            qrg(:, nx)   = qrg(:, 0)

            Nrg(:, 0)    = Nrg(:, 0) + Nrg(:, nx)
            Nrg(:, nx-1) = Nrg(:, nx-1) + Nrg(:, -1)
            Nrg(:, -1)   = Nrg(:, nx-1)
            Nrg(:, nx)   = Nrg(:, 0)

            ! apply free slip boundary condition
            prec_volg(0,  :) = two * prec_volg(0,  :)
            prec_volg(nz, :) = two * prec_volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            prec_volg(1,    :) = prec_volg(1,    :) + prec_volg(-1,   :)
            prec_volg(nz-1, :) = prec_volg(nz-1, :) + prec_volg(nz+1, :)

            prec_tbuoyg(0,  :) = two * prec_tbuoyg(0,  :)
            prec_tbuoyg(nz, :) = two * prec_tbuoyg(nz, :)
            prec_tbuoyg(1,    :) = prec_tbuoyg(1,    :) + prec_tbuoyg(-1,   :)
            prec_tbuoyg(nz-1, :) = prec_tbuoyg(nz-1, :) + prec_tbuoyg(nz+1, :)
            prec_tbuoyg(0:nz, :) = prec_tbuoyg(0:nz, :) / volg(0:nz, :) ! Note to divide by volg, not prec_volg!
            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            prec_tbuoyg(-1,   :) = two * prec_tbuoyg(0,  :) - prec_tbuoyg(1, :)
            prec_tbuoyg(nz+1, :) = two * prec_tbuoyg(nz, :) - prec_tbuoyg(nz-1, :)

            qrg(0,  :) = two * qrg(0,  :)
            qrg(nz, :) = two * qrg(nz, :)
            qrg(1,    :) = qrg(1,    :) + qrg(-1,   :)
            qrg(nz-1, :) = qrg(nz-1, :) + qrg(nz+1, :)
            qrg(0:nz, :) = qrg(0:nz, :) / volg(0:nz, :) ! Note to divide by volg, not prec_volg!
            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            qrg(-1,   :) = two * qrg(0,  :) - qrg(1, :)
            qrg(nz+1, :) = two * qrg(nz, :) - qrg(nz-1, :)

            Nrg(0,  :) = two * Nrg(0,  :)
            Nrg(nz, :) = two * Nrg(nz, :)
            Nrg(1,    :) = Nrg(1,    :) + Nrg(-1,   :)
            Nrg(nz-1, :) = Nrg(nz-1, :) + Nrg(nz+1, :)
            Nrg(0:nz, :) = Nrg(0:nz, :) / volg(0:nz, :) ! Note to divide by volg, not prec_volg!
            ! extrapolate to halo grid points (needed to compute
            ! z derivative used for the time step)
            Nrg(-1,   :) = two * Nrg(0,  :) - Nrg(1, :)
            Nrg(nz+1, :) = two * Nrg(nz, :) - Nrg(nz-1, :)

            ! sum halo contribution into internal cells
            ! (be aware that halo cell contribution at upper boundary
            ! are added to cell nz)
            prec_nparg(0,    :) = prec_nparg(0,    :) + prec_nparg(-1, :)
            prec_nparg(nz-1, :) = prec_nparg(nz-1, :) + prec_nparg(nz, :)

            ! sanity check
            if (sum(prec_nparg(0:nz-1, :)) /= n_prec_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

            call stop_timer(prec_par2grid_timer)

        end subroutine prec_par2grid


        ! Interpolate the gridded quantities to the parcels
        ! @param[inout] vel is the parcel velocity
        ! @param[in] add contributions, i.e. do not reset parcel quantities to zero before doing grid2par.
        !            (optional)
        subroutine prec_grid2par(vel, add)
            double precision,     intent(inout) :: vel(:, :)
            logical, optional, intent(in)       :: add
            double precision                    :: points(2)
            integer                             :: n, p, l

            call start_timer(prec_grid2par_timer)

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    !$omp parallel default(shared)
                    !$omp do private(n)
                    do n = 1, n_prec_parcels
                        vel(:, n) = zero
                    enddo
                    !$omp end do
                    !$omp end parallel
               endif
            else
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_prec_parcels
                    vel(:, n) = zero
                enddo
                !$omp end do
                !$omp end parallel
            endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, is, js, weights)
            do n = 1, n_prec_parcels
                points = prec_parcels%position(:, n)

                call apply_periodic_bc(points)

                ! get interpolation weights and mesh indices
                call bilinear(points(:), is, js, weights)

                do l = 1,2
                    vel(l, n) = vel(l, n) &
                              + sum(weights * velog(js:js+1, is:is+1, l))
                end do
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(prec_grid2par_timer)

        end subroutine prec_grid2par


        ! Interpolate the gridded quantities to the parcels without resetting
        ! their values to zero before doing grid2par.
        ! @param[inout] vel is the parcel velocity
        subroutine prec_grid2par_add(vel)
            double precision,       intent(inout) :: vel(:, :)

            call prec_grid2par(vel, add=.true.)

        end subroutine prec_grid2par_add

end module prec_parcel_interpl
