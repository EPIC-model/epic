! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12, thousand
    use merge_sort
    use parameters, only : extent, lower, vcell, vmin, nx, nz, vdomaini
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid
    use omp_lib
    use physics, only : peref, ape_calculation
    use timer, only : start_timer, stop_timer
    use ape_density, only : ape_den
    implicit none

    integer :: parcel_stats_timer

#ifndef NDEBUG
    double precision, parameter :: thres = thousand * epsilon(zero)
#endif

    ! ape   : available potential energy
    ! ke    : kinetic energy
    double precision :: ape, ke

    ! en : enstrophy
    double precision :: en

    integer :: n_small

    ! avg_lam : mean aspect ratio over all parcels
    ! avg_vol : mean volume over all parcels
    ! std_lam : standard deviation of aspect ratio
    ! std_vol : standard deviation of volume
    double precision :: avg_lam, avg_vol
    double precision :: std_lam, std_vol

    double precision :: sum_vol

    ! rms vorticity
    double precision :: rms_zeta(3)

    ! buoyancy extrema
    double precision :: bmin, bmax

    contains

        ! Compute the reference potential energy
        subroutine calculate_peref
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            if (.not. ape_calculation == "sorting") then
                peref = zero
                return
            endif

            call start_timer(parcel_stats_timer)

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = f12 / (extent(1) * extent(2))
            zmean = gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * (parcels%volume(ii(n-1)) + parcels%volume(ii(n)))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo

            ! divide by domain volume to get domain-averaged peref
            peref = peref * vdomaini

            call stop_timer(parcel_stats_timer)
        end subroutine calculate_peref


        ! Calculate all parcel related diagnostics
        subroutine calculate_parcel_diagnostics(velocity)
            double precision, intent(in) :: velocity(:, :)
            integer          :: n
            double precision :: b, z, vel(3), vol, zmin, vor(3), pe
            double precision :: evals(3), lam, lsum, l2sum, v2sum

            call start_timer(parcel_stats_timer)

            ! reset
            ke = zero
            ape = zero
            en = zero

            lsum = zero
            l2sum = zero
            v2sum = zero

            rms_zeta = zero

            n_small = zero

            zmin = lower(3)

            avg_lam = zero
            avg_vol = zero
            std_lam = zero
            std_vol = zero
            sum_vol = zero
            bmin = zero
            bmax = zero

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z, evals, lam, vor) &
            !$omp& reduction(+: ke, ape, lsum, l2sum, sum_vol, v2sum, n_small, rms_zeta, en) &
            !$omp& reduction(-: pe) reduction(max: bmax) reduction(min: bmin)
            do n = 1, n_parcels

                vel = velocity(:, n)
                vor = parcels%vorticity(:, n)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(3, n)

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2 + vel(3) ** 2) * vol

                if (ape_calculation == 'sorting') then
                    ! potential energy using sorting approach
                    pe = pe - b * (z - zmin) * vol
                else if (ape_calculation == 'ape density') then
                    ape = ape + ape_den(b, z) * vol
                endif

                ! enstrophy
                en = en + (vor(1) ** 2 + vor(2) ** 2 + vor(3) ** 2) * vol

                evals = get_eigenvalues(parcels%B(:, n), parcels%volume(n))
                lam = get_aspect_ratio(evals)

                lsum = lsum + lam
                l2sum = l2sum + lam ** 2

                sum_vol = sum_vol + vol
                v2sum = v2sum + vol ** 2

                if (vol <= vmin) then
                    n_small = n_small + 1
                endif

                if (bmax < b) then
                    bmax = b
                endif

                if (bmin > b) then
                    bmin = b
                endif

#ifndef NDEBUG
                !$omp critical
                if (abs(get_determinant(parcels%B(:, n), vol) / get_abc(vol) ** 2 - one) > thres) then
                    print *, "Parcel determinant not preserved!"
                    stop
                endif
                !$omp end critical
#endif
                rms_zeta = rms_zeta + vol * vor ** 2

            enddo
            !$omp end do
            !$omp end parallel

            ! divide by domain volume to get domain-averaged quantities
            ke = f12 * ke * vdomaini

            if (ape_calculation == 'sorting') then
                ape = pe * vdomaini - peref
            else if (ape_calculation == 'ape density') then
                ape = ape * vdomaini
            endif

            en = f12 * en * vdomaini

            avg_lam = lsum / dble(n_parcels)
            std_lam = dsqrt(abs(l2sum / dble(n_parcels) - avg_lam ** 2))

            rms_zeta = dsqrt(rms_zeta / sum_vol)

            avg_vol = sum_vol / dble(n_parcels)
            std_vol = dsqrt(abs(v2sum / dble(n_parcels) - avg_vol ** 2))

            call stop_timer(parcel_stats_timer)

        end subroutine calculate_parcel_diagnostics

end module parcel_diagnostics
