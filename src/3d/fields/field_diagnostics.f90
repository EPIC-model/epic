! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli
    use constants, only : f12
    use fields
    use timer, only : start_timer, stop_timer
    use mpi_layout, only : box
    use mpi_communicator
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    integer :: field_stats_timer

    ! Array indices of field stats array
    integer, parameter :: IDX_RMS_V         = 1, &  ! rms volume error
                          IDX_AVG_NPAR      = 2, &  ! average num parcels per cell
                          IDX_AVG_NSPAR     = 3, &  ! average num small parcels per cell
                          IDX_KEG           = 4, &  ! kinetic energy calculated on the grid
                          IDX_ABSERR_V      = 5, &  ! max absolute normalised volume error
                          IDX_MAX_NPAR      = 6, &  ! max num parcels per cell
                          IDX_MIN_NPAR      = 7     ! min num parcels per cell

    double precision :: field_stats(IDX_MIN_NPAR)

    contains

        ! Note: Only the MPI root has the valid data after
        ! this operation.
        subroutine calculate_field_diagnostics
            integer          :: lo(3), hi(3)

            call start_timer(field_stats_timer)

            lo = box%lo
            hi = box%hi

            !
            ! calculate locally
            !

            ! do not take halo cells into account
            field_stats(IDX_RMS_V) = sum((volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)) - vcell) ** 2)

            field_stats(IDX_ABSERR_V) = maxval(abs(volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))  - vcell)) * vcelli

            field_stats(IDX_MAX_NPAR) = maxval(nparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1)))

            field_stats(IDX_MIN_NPAR) = minval(nparg(lo(3):hi(3)-1, lo(2):hi(2), lo(1):hi(1)))

            field_stats(IDX_AVG_NPAR) = sum(nparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))) * ncelli

            field_stats(IDX_AVG_NSPAR) = sum(nsparg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1))) * ncelli

            field_stats(IDX_KEG) = f12 * sum( volg(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)) *         &
                                            (velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 1) ** 2   &
                                           + velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 2) ** 2   &
                                           + velog(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1), 3) ** 2))


            !
            ! do communication
            !
            call mpi_blocking_reduce(field_stats(IDX_RMS_V:IDX_KEG), MPI_SUM)

            call mpi_blocking_reduce(field_stats(IDX_ABSERR_V:IDX_MAX_NPAR), MPI_MAX)

            call mpi_blocking_reduce(field_stats(IDX_MIN_NPAR), MPI_MIN)

            !
            ! final calculations
            !
            field_stats(IDX_RMS_V) = dsqrt(field_stats(IDX_RMS_V) * ngridi) * vcelli

            call stop_timer(field_stats_timer)

        end subroutine calculate_field_diagnostics

end module field_diagnostics
