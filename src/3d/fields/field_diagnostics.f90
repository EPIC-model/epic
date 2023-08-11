! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli, vdomaini
    use constants, only : f12, f14
    use fields
    use mpi_timer, only : start_timer, stop_timer
    use mpi_layout, only : box
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use physics, only : ape_calculation
    use ape_density, only : ape_den
    implicit none

    integer :: field_stats_timer

    ! Array indices of field stats array
    integer, parameter :: IDX_RMS_V         = 1,  & ! rms volume error
                          IDX_AVG_NPAR      = 2,  & ! average num parcels per cell
                          IDX_AVG_NSPAR     = 3,  & ! average num small parcels per cell
                          IDX_KEG           = 4,  & ! domain-averaged kinetic energy calculated on the grid
                          IDX_APEG          = 5,  & ! domain-averaged available potential energy on the grid
                          IDX_ENG           = 6,  & ! domain-averaged enstrophy calculated on the grid
                          IDX_ABSERR_V      = 7,  & ! max absolute normalised volume error
                          IDX_MAX_NPAR      = 8,  & ! max num parcels per cell
                          IDX_MAX_BUOY      = 9,  & ! max gridded buoyancy
                          IDX_MIN_NPAR      = 10, & ! min num parcels per cell
                          IDX_MIN_BUOY      = 11    ! min gridded buoyancy

    double precision :: field_stats(IDX_MIN_BUOY)

    contains

        ! Note: Only the MPI root has the valid data after
        ! this operation.
        subroutine calculate_field_diagnostics
            integer          :: lo(3), hi(3)
            double precision :: z(0:nz)
            integer          :: ix, iy, iz

            call start_timer(field_stats_timer)

            lo = box%lo
            hi = box%hi

            !
            ! calculate locally
            !

            ! do not take halo cells into account
            field_stats(IDX_RMS_V) = sum((volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)) - vcell) ** 2)

            field_stats(IDX_ABSERR_V) = maxval(abs(volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))  - vcell)) * vcelli

            ! vertically only include 0 (= lo(3) to nz-1 (= hi(3) - 1)
            field_stats(IDX_MAX_NPAR) = maxval(nparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1)))

            field_stats(IDX_MIN_NPAR) = minval(nparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1)))

            field_stats(IDX_AVG_NPAR) = sum(nparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1))) * ncelli

            field_stats(IDX_AVG_NSPAR) = sum(nsparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1))) * ncelli

            ! do not take halo cells into account
            field_stats(IDX_MIN_BUOY) = minval(tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)))
            field_stats(IDX_MAX_BUOY) = maxval(tbuoyg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)))

            ! use half weights for boundary grid points
            field_stats(IDX_KEG) = f12 * sum( volg(1:nz-1, lo(2):hi(2), lo(1):hi(1))           &
                                          * (velog(1:nz-1, lo(2):hi(2), lo(1):hi(1), 1) ** 2   &
                                           + velog(1:nz-1, lo(2):hi(2), lo(1):hi(1), 2) ** 2   &
                                           + velog(1:nz-1, lo(2):hi(2), lo(1):hi(1), 3) ** 2)) &
                                 + f14 * sum(volg( 0,  lo(2):hi(2), lo(1):hi(1))               &
                                          * (velog(0,  lo(2):hi(2), lo(1):hi(1), 1) ** 2       &
                                           + velog(0,  lo(2):hi(2), lo(1):hi(1), 2) ** 2       &
                                           + velog(0,  lo(2):hi(2), lo(1):hi(1), 3) ** 2))     &
                                 + f14 * sum( volg(nz, lo(2):hi(2), lo(1):hi(1))               &
                                          * (velog(nz, lo(2):hi(2), lo(1):hi(1), 1) ** 2       &
                                           + velog(nz, lo(2):hi(2), lo(1):hi(1), 2) ** 2       &
                                           + velog(nz, lo(2):hi(2), lo(1):hi(1), 3) ** 2))

            ! divide by domain volume to get domain-averaged kinetic energy
            field_stats(IDX_KEG) = field_stats(IDX_KEG) * vdomaini

            field_stats(IDX_ENG) = f12 * sum( volg(1:nz-1, lo(2):hi(2), lo(1):hi(1))           &
                                          * (vortg(1:nz-1, lo(2):hi(2), lo(1):hi(1), 1) ** 2   &
                                           + vortg(1:nz-1, lo(2):hi(2), lo(1):hi(1), 2) ** 2   &
                                           + vortg(1:nz-1, lo(2):hi(2), lo(1):hi(1), 3) ** 2)) &
                                 + f14 * sum(volg( 0,  lo(2):hi(2), lo(1):hi(1))               &
                                          * (vortg(0,  lo(2):hi(2), lo(1):hi(1), 1) ** 2       &
                                           + vortg(0,  lo(2):hi(2), lo(1):hi(1), 2) ** 2       &
                                           + vortg(0,  lo(2):hi(2), lo(1):hi(1), 3) ** 2))     &
                                 + f14 * sum( volg(nz, lo(2):hi(2), lo(1):hi(1))               &
                                          * (vortg(nz, lo(2):hi(2), lo(1):hi(1), 1) ** 2       &
                                           + vortg(nz, lo(2):hi(2), lo(1):hi(1), 2) ** 2       &
                                           + vortg(nz, lo(2):hi(2), lo(1):hi(1), 3) ** 2))

            ! divide by domain volume to get domain-averaged enstrophy
            field_stats(IDX_ENG) = field_stats(IDX_ENG) * vdomaini

            if (ape_calculation == 'ape density') then
                do iz = 0, nz
                    z(iz) = lower(3) + dble(iz) * dx(3)
                enddo

                field_stats(IDX_APEG) = zero
                do ix = lo(1), hi(1)
                    do iy = lo(2), hi(2)
                        field_stats(IDX_APEG) = field_stats(IDX_APEG) &
                                      +   sum(volg(1:nz-1, iy, ix) * ape_den(tbuoyg(1:nz-1, iy, ix), z(1:nz-1))) &
                                      + f12 * volg(0,      iy, ix) * ape_den(tbuoyg(0,      iy, ix), z(0))       &
                                      + f12 * volg(nz,     iy, ix) * ape_den(tbuoyg(nz,     iy, ix), z(nz))
                    enddo
                enddo

                field_stats(IDX_APEG) = field_stats(IDX_APEG) * vdomaini
            endif

            !
            ! do communication
            !
            call mpi_blocking_reduce(field_stats(IDX_RMS_V:IDX_ENG), MPI_SUM, world)

            call mpi_blocking_reduce(field_stats(IDX_ABSERR_V:IDX_MAX_BUOY), MPI_MAX, world)

            call mpi_blocking_reduce(field_stats(IDX_MIN_NPAR:IDX_MIN_BUOY), MPI_MIN, world)

            !
            ! final calculations
            !
            field_stats(IDX_RMS_V) = dsqrt(field_stats(IDX_RMS_V) * ngridi) * vcelli

            call stop_timer(field_stats_timer)

        end subroutine calculate_field_diagnostics

end module field_diagnostics
