! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli
    use constants, only : f12, f14
    use fields
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_stats_timer

    double precision :: rms_v,      &       ! rms volume error
                        abserr_v,   &       ! max absolute normalised volume error
                        max_npar,   &       ! max num parcels per cell
                        min_npar,   &       ! min num parcels per cell
                        avg_npar,   &       ! average num parcels per cell
                        avg_nspar,  &       ! average num small parcels per cell
                        keg                 ! kinetic energy calculated on the grid
    contains

        subroutine calculate_field_diagnostics
            double precision :: sqerrsum

            call start_timer(field_stats_timer)

            ! do not take halo cells into account
            sqerrsum = sum((volg(0:nz, :, :) - vcell) ** 2)
            rms_v = dsqrt(sqerrsum * ngridi) * vcelli

            abserr_v = maxval(abs(volg(0:nz, :, :)  - vcell)) * vcelli

            max_npar = maxval(nparg(0:nz-1, :, :))

            min_npar = minval(nparg(0:nz-1, :, :))

            avg_npar = sum(nparg(0:nz-1, :, :)) * ncelli

            avg_nspar = sum(nsparg(0:nz-1, :, :)) * ncelli

            ! use half weights for boundary grid points
            keg = f12 * sum(volg(1:nz-1, :, :) * ( velog(1:nz-1, :, :, 1) ** 2   &
                                                 + velog(1:nz-1, :, :, 2) ** 2   &
                                                 + velog(1:nz-1, :, :, 3) ** 2)) &
                + f14 * sum(volg(0,  :, :) * ( velog(0,  :, :, 1) ** 2           &
                                             + velog(0,  :, :, 2) ** 2           &
                                             + velog(0,  :, :, 3) ** 2))         &
                + f14 * sum(volg(nz, :, :) * ( velog(nz, :, :, 1) ** 2           &
                                             + velog(nz, :, :, 2) ** 2           &
                                             + velog(nz, :, :, 3) ** 2))

            call stop_timer(field_stats_timer)

        end subroutine calculate_field_diagnostics

end module field_diagnostics
