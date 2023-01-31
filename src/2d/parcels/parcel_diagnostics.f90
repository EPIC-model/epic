! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower, vcell, vmin, nx, nz, vdomaini
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse
    use omp_lib
    use physics, only : peref, ape_calculation
    use timer, only : start_timer, stop_timer
    use ape_density, only : ape_den
    implicit none

    integer :: parcel_stats_timer

    ! pe    : domain-averaged potential energy
    ! ape   : domain-average available potential energy
    ! ke    : domain-averaged kinetic energy
    double precision :: pe, ape, ke

    integer :: n_small

    ! avg_lam : mean aspect ratio over all parcels
    ! avg_vol : mean volume over all parcels
    ! std_lam : standard deviation of aspect ratio
    ! std_vol : standard deviation of volume
    double precision :: avg_lam, avg_vol
    double precision :: std_lam, std_vol

    ! rms vorticity
    double precision :: rms_zeta

    ! min and max buoyancy
    double precision :: bmin, bmax

    ! min and max vorticity
    double precision :: vormin, vormax

#ifdef ENABLE_DIAGNOSE
    ! buoyancy weighted first and second moments
    double precision :: xb_bar, x2b_bar
    double precision :: zb_bar, z2b_bar
    double precision :: xzb_bar

    ! vorticity weighted first and second moments
    double precision :: xv_bar, x2v_bar
    double precision :: zv_bar, z2v_bar
    double precision :: xzv_bar
#endif

    contains

        ! Compute the reference potential energy
        subroutine calculate_peref
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            call start_timer(parcel_stats_timer)

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = f12 / extent(1)
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
            double precision :: velocity(:, :)
            integer          :: n
            double precision :: b, z, vel(2), vol, zmin
            double precision :: eval, lam, B22, lsum, l2sum, vsum, v2sum

            call start_timer(parcel_stats_timer)

            ! reset
            ke = zero
            ape = zero
            pe = zero

            ! find extrema outside OpenMP loop, we can integrate it later;
            ! this way the result is reproducible
            bmin = minval(parcels%buoyancy(1:n_parcels))
            bmax = maxval(parcels%buoyancy(1:n_parcels))
            vormin = minval(parcels%vorticity(1:n_parcels))
            vormax = maxval(parcels%vorticity(1:n_parcels))

            lsum = zero
            l2sum = zero
            vsum = zero
            v2sum = zero

            rms_zeta = zero

            n_small = zero

            zmin = lower(2)

            avg_lam = zero
            avg_vol = zero
            std_lam = zero
            std_vol = zero

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z, eval, lam, B22) &
            !$omp& reduction(+: ke, ape, lsum, l2sum, vsum, v2sum, n_small, rms_zeta) &
            !$omp& reduction(-: pe)
            do n = 1, n_parcels

                vel = velocity(:, n)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(2, n)

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                if (ape_calculation == 'sorting') then
                    ! potential energy using sorting approach
                    pe = pe - b * (z - zmin) * vol
                else if (ape_calculation == 'ape density') then
                    ape = ape + ape_den(b, z) * vol
                endif

                B22 = get_B22(parcels%B(1, n), parcels%B(2, n), vol)
                eval = get_eigenvalue(parcels%B(1, n), parcels%B(2, n), B22)
                lam = get_aspect_ratio(eval, vol)

                lsum = lsum + lam
                l2sum = l2sum + lam ** 2

                vsum = vsum + vol
                v2sum = v2sum + vol ** 2

                if (vol <= vmin) then
                    n_small = n_small + 1
                endif

                rms_zeta = rms_zeta + vol * parcels%vorticity(n) ** 2

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

            avg_lam = lsum / dble(n_parcels)
            std_lam = dsqrt(abs(l2sum / dble(n_parcels) - avg_lam ** 2))

            rms_zeta = dsqrt(rms_zeta / vsum)

            avg_vol = vsum / dble(n_parcels)
            std_vol = dsqrt(abs(v2sum / dble(n_parcels) - avg_vol ** 2))


#ifdef ENABLE_DIAGNOSE
            call straka_diagnostics
#endif

            call stop_timer(parcel_stats_timer)

        end subroutine calculate_parcel_diagnostics


#ifdef ENABLE_DIAGNOSE
        ! Straka density current test case diagnostics
        subroutine straka_diagnostics
            integer          :: n
            double precision :: xbv, x2bv, zbv, z2bv, xzbv
            double precision :: xvv, x2vv, zvv, z2vv, xzvv
            double precision :: bvsum, vvsum, bv, vv

            ! reset
            xb_bar = zero
            zb_bar = zero
            x2b_bar = zero
            z2b_bar = zero
            xzb_bar = zero

            xv_bar = zero
            zv_bar = zero
            x2v_bar = zero
            z2v_bar = zero
            xzv_bar = zero


            xbv = zero
            x2bv = zero
            zbv = zero
            z2bv = zero
            xzbv = zero
            xvv = zero
            x2vv = zero
            zvv = zero
            z2vv = zero
            xzvv = zero

            bvsum = zero
            vvsum = zero
            bv = zero
            vv = zero

            !$omp parallel default(shared)
            !$omp do private(n, bv, vv) &
            !$omp& reduction(+: vvsum, bvsum, xbv, zbv, x2bv, z2bv, xzbv, xvv, zvv, x2vv, z2vv, xzvv)
            do n = 1, n_parcels
                ! we only use the upper half in horizontal direction
                if (parcels%position(1, n) >= 0) then
                    bv = parcels%buoyancy(n) * parcels%volume(n)
                    bvsum = bvsum + bv
                    xbv = xbv + bv * parcels%position(1, n)
                    zbv = zbv + bv * parcels%position(2, n)

                    x2bv = x2bv + bv * parcels%position(1, n) ** 2
                    z2bv = z2bv + bv * parcels%position(2, n) ** 2
                    xzbv = xzbv + bv * parcels%position(1, n) * parcels%position(2, n)


                    vv = parcels%vorticity(n) * parcels%volume(n)
                    vvsum = vvsum + vv
                    xvv = xvv + vv * parcels%position(1, n)
                    zvv = zvv + vv * parcels%position(2, n)

                    x2vv = x2vv + vv * parcels%position(1, n) ** 2
                    z2vv = z2vv + vv * parcels%position(2, n) ** 2
                    xzvv = xzvv + vv * parcels%position(1, n) * parcels%position(2, n)
                endif
            enddo
            !$omp end do
            !$omp end parallel

            ! we do not need to divide by the number of of involved
            ! parcels since whe divide by the sums "bvsum" or "vvsum"
            ! what should also be averages (i.e. divided by the number of
            ! involved parcels)

            ! make sure we do not divide by zero
            if (dabs(bvsum) < epsilon(zero)) then
                bvsum = epsilon(zero)
            else
                bvsum = one / bvsum
            endif

            if (dabs(vvsum) < epsilon(zero)) then
                vvsum = epsilon(zero)
            else
                vvsum = one / vvsum
            endif

            xb_bar = xbv * bvsum
            zb_bar = zbv * bvsum

            x2b_bar = x2bv * bvsum
            z2b_bar = z2bv * bvsum

            xzb_bar = xzbv * bvsum - xb_bar * zb_bar

            xv_bar = xvv * vvsum
            zv_bar = zvv * vvsum

            x2v_bar = x2vv * vvsum
            z2v_bar = z2vv * vvsum

            xzv_bar = xzvv * vvsum - xv_bar * zv_bar
        end subroutine straka_diagnostics
#endif
end module parcel_diagnostics
