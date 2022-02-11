! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli
    use fields
    implicit none

    double precision :: rms_v,      &       ! rms volume error
                        abserr_v,   &       ! max absolute normalised volume error
                        max_npar,   &       ! max num parcels per cell
                        min_npar,   &       ! min num parcels per cell
                        avg_npar,   &       ! average num parcels per cell
                        avg_nspar           ! average num small parcels per cell
    contains

        subroutine calculate_field_diagnostics
            double precision :: sqerrsum

            ! do not take halo cells into account
            sqerrsum = sum((volg(0:nz, :, :) - vcell) ** 2)
            rms_v = dsqrt(sqerrsum * ngridi) * vcelli

            abserr_v = maxval(abs(volg(0:nz, :, :)  - vcell)) * vcelli

            max_npar = maxval(nparg(0:nz-1, :, :))

            min_npar = minval(nparg(0:nz-1, :, :))

            avg_npar = sum(nparg(0:nz-1, :, :)) * ncelli

            avg_nspar = sum(nsparg(0:nz-1, :, :)) * ncelli

        end subroutine calculate_field_diagnostics

end module field_diagnostics
