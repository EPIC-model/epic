! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module surface_parcel_diagnostics
    use constants, only : zero, one, f12
    use parameters, only : amin, nx, nz, asurfi
    use surface_parcel_container, only : surface_parcel_container_type
    use parcel_ellipse
    use omp_lib
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: surf_parcel_stats_timer

    ! ke    : domain-averaged kinetic energy
    double precision :: ke(2)

    integer :: n_small(2)

    ! avg_lam  : mean aspect ratio over all parcels
    ! avg_area : mean area over all parcels
    ! std_lam  : standard deviation of aspect ratio
    ! std_area : standard deviation of area
    double precision :: avg_lam(2), avg_area(2)
    double precision :: std_lam(2), std_area(2)

    ! rms vorticity
    double precision :: rms_xi(2)
    double precision :: rms_eta(2)
    double precision :: rms_zeta(2)

    ! min and max buoyancy
    double precision :: bmin(2), bmax(2)

    contains

        ! Calculate all parcel related diagnostics
        subroutine calculate_surface_parcel_diagnostics(s_parcels, n_par, which, velocity)
            type(surface_parcel_container_type), intent(in) :: s_parcels
            integer,                             intent(in) :: n_par
            double precision,                    intent(in) :: velocity(:, :)
            character(2),                        intent(in) :: which
            integer                                         :: n, j
            double precision                                :: vel(2), area
            double precision                                :: eval, lam, lsum
            double precision                                :: l2sum, asum, a2sum

            call start_timer(surf_parcel_stats_timer)

            j = 1
            if (which == 'up') then
                j = 2
            endif

            ! reset
            ke(j) = zero

            ! find extrema outside OpenMP loop, we can integrate it later;
            ! this way the result is reproducible
            bmin(j) = minval(s_parcels%buoyancy(1:n_par))
            bmax(j) = maxval(s_parcels%buoyancy(1:n_par))

            lsum = zero
            l2sum = zero
            asum = zero
            a2sum = zero

            rms_xi(j) = zero
            rms_eta(j) = zero
            rms_zeta(j) = zero

            n_small(j) = zero

            avg_lam(j) = zero
            avg_area(j) = zero
            std_lam(j) = zero
            std_area(j) = zero

            !$omp parallel default(shared)
            !$omp do private(n, vel, area, eval, lam) &
            !$omp& reduction(+: ke, lsum, l2sum, asum, a2sum, n_small, rms_xi, rms_eta, rms_zeta)
            do n = 1, n_par

                vel  = velocity(:, n)
                area = s_parcels%area(n)

                ! kinetic energy
                ke(j) = ke(j) + (vel(1) ** 2 + vel(2) ** 2) * area

                eval = get_eigenvalue(s_parcels%B(:, n))
                lam = get_aspect_ratio(eval, area)

                lsum = lsum + lam
                l2sum = l2sum + lam ** 2

                asum = asum + area
                a2sum = a2sum + area ** 2

                if (area <= amin) then
                    n_small(j) = n_small(j) + 1
                endif


                rms_xi(j)   = rms_xi(j)   + area * s_parcels%vorticity(1, n) ** 2
                rms_eta(j)  = rms_eta(j)  + area * s_parcels%vorticity(2, n) ** 2
                rms_zeta(j) = rms_zeta(j) + area * s_parcels%vorticity(3, n) ** 2

            enddo
            !$omp end do
            !$omp end parallel

            ! divide by domain volume to get domain-averaged quantities
            ke(j) = f12 * ke(j) * asurfi

            avg_lam(j) = lsum / dble(n_par)
            std_lam(j) = dsqrt(abs(l2sum / dble(n_par) - avg_lam(j) ** 2))

            rms_xi(j)   = dsqrt(rms_xi(j)   / asum)
            rms_eta(j)  = dsqrt(rms_eta(j)  / asum)
            rms_zeta(j) = dsqrt(rms_zeta(j) / asum)

            avg_area(j) = asum / dble(n_par)
            std_area(j) = dsqrt(abs(a2sum / dble(n_par) - avg_area(j) ** 2))


            call stop_timer(surf_parcel_stats_timer)

        end subroutine calculate_surface_parcel_diagnostics

end module surface_parcel_diagnostics
