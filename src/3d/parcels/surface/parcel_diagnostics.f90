! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower, acell, amin, nx, ny
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse
    use omp_lib
    use physics, only : csq
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: parcel_stats_timer

    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: pe, ke

    integer :: n_small

    ! avg_lam  : mean aspect ratio over all parcels
    ! avg_area : mean area over all parcels
    ! std_lam  : standard deviation of aspect ratio
    ! std_area : standard deviation of area
    double precision :: avg_lam, avg_area
    double precision :: std_lam, std_area

    ! rms potential vorticity
    double precision :: rms_q

    ! min and max potential vorticity
    double precision :: min_q, max_q

    contains

        ! Calculate all parcel related diagnostics
        subroutine calculate_parcel_diagnostics(velocity)
            double precision :: velocity(:, :)
            integer          :: n
            double precision :: vel(2), area, mass
            double precision :: eval, lam, lsum, l2sum, asum, a2sum

            call start_timer(parcel_stats_timer)

            ! reset
            ke = zero
            pe = zero

            ! find extrema outside OpenMP loop, we can integrate it later;
            ! this way the result is reproducible
            min_q = minval(parcels%q(1:n_parcels))
            max_q = maxval(parcels%q(1:n_parcels))

            lsum = zero
            l2sum = zero
            asum = zero
            a2sum = zero

            rms_q = zero

            n_small = zero

            avg_lam = zero
            avg_area = zero
            std_lam = zero
            std_area = zero

            !$omp parallel default(shared)
            !$omp do private(n, vel, area, mass, eval, lam) &
            !$omp& reduction(+: ke, pe, lsum, l2sum, asum, a2sum, n_small, rms_q)
            do n = 1, n_parcels

                vel = velocity(:, n)
                area = parcels%area(n)
                mass = parcels%mass(n)

                ! kinetic energy divided by mean depth H:
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * area * (mass + one)

                ! potential energy divided by mean depth H:
                pe = pe + area * csq * mass ** 2

                eval = get_eigenvalue(parcels%B(:, n))
                lam = get_aspect_ratio(eval, area)

                lsum = lsum + lam
                l2sum = l2sum + lam ** 2

                asum = asum + area
                a2sum = a2sum + area ** 2

                if (area <= amin) then
                    n_small = n_small + 1
                endif

                rms_q = rms_q + area * parcels%q(n) ** 2

            enddo
            !$omp end do
            !$omp end parallel

            ke = f12 * ke
            pe = f12 * pe

            avg_lam = lsum / dble(n_parcels)
            std_lam = dsqrt(abs(l2sum / dble(n_parcels) - avg_lam ** 2))

            rms_q = dsqrt(rms_q / asum)

            avg_area = asum / dble(n_parcels)
            std_area = dsqrt(abs(a2sum / dble(n_parcels) - avg_area ** 2))

            call stop_timer(parcel_stats_timer)

        end subroutine calculate_parcel_diagnostics

end module parcel_diagnostics
