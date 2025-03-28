! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use datatypes, only : int64
    use constants, only : zero, one, f12, thousand
    use parameters, only : extent, lower, vcell, vmin, nx, nz, vdomaini
    use parcels_mod, only : parcels
    use parcel_split_mod, only : n_parcel_splits
    use parcel_merging, only : n_parcel_merges, n_big_close, n_way_parcel_mergers
    use omp_lib
    use physics, only : ape_calculation
    use ape_density, only : ape_den
    use mpi_timer, only : start_timer, stop_timer
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    integer :: parcel_stats_timer

#ifndef NDEBUG
    double precision, parameter :: thres = thousand * epsilon(zero)
#endif

    integer, parameter :: IDX_APE       =  1, & ! available potential energy
                          IDX_KE        =  2, & ! kinetic energy
                          IDX_N_SMALL   =  3, & ! number of small parcels (V_i < V_min)
                          IDX_AVG_LAM   =  4, & ! mean aspect ratio over all parcels
                          IDX_AVG_VOL   =  5, & ! mean volume over all parcels
                          IDX_STD_LAM   =  6, & ! standard deviation of aspect ratio
                          IDX_STD_VOL   =  7, & ! standard deviation of volume
                          IDX_SUM_VOL   =  8, & ! total volume
                          IDX_RMS_XI    =  9, & ! rms value of x-vorticity
                          IDX_RMS_ETA   = 10, & ! rms value of y-vorticity
                          IDX_RMS_ZETA  = 11, & ! rms value of z-vorticity
                          IDX_NTOT_PAR  = 12, & ! total number of parcels
                          IDX_ENSTROPHY = 13, & ! enstrophy
                          IDX_NSPLITS   = 14, & ! number of parcel splits since last write
                          IDX_NBIG_ICLO = 15, & ! number of big iclo neighbours (merging)
                          IDX_NMERGES   = 16, & ! number of parcel merges since last write
                          IDX_MIN_BUOY  = 17, & ! minimum parcel buoyancy
                          IDX_MAX_BUOY  = 18    ! maximum parcel buoyancy

    double precision     :: parcel_stats(IDX_MAX_BUOY)
    integer(kind=int64)  :: parcel_merge_stats(size(n_way_parcel_mergers))

    contains

        ! Calculate all parcel related diagnostics
        subroutine calculate_parcel_diagnostics
            integer          :: n
            double precision :: b, z, vel(3), vol, vor(3)
            double precision :: ntoti, bmin, bmax
            double precision :: evals(3), lam

            call start_timer(parcel_stats_timer)

            ! reset
            parcel_stats = zero
            bmin = zero
            bmax = zero

            parcel_stats(IDX_NTOT_PAR) = dble(parcels%local_num)

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z, evals, lam, vor) &
            !$omp& reduction(+: parcel_stats) &
            !$omp& reduction(max: bmax) reduction(min: bmin)
            do n = 1, parcels%local_num

                vel = parcels%delta_pos(:, n)
                vor = parcels%vorticity(:, n)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(3, n)

                ! kinetic energy
                parcel_stats(IDX_KE) = parcel_stats(IDX_KE) + (vel(1) ** 2 + vel(2) ** 2 + vel(3) ** 2) * vol

                if (ape_calculation == 'ape density') then
                    parcel_stats(IDX_APE) = parcel_stats(IDX_APE) + ape_den(b, z) * vol
                endif

                ! enstrophy
                parcel_stats(IDX_ENSTROPHY) = parcel_stats(IDX_ENSTROPHY) &
                                            + (vor(1) ** 2 + vor(2) ** 2 + vor(3) ** 2) * vol

                evals = parcels%get_eigenvalues(n)
                lam = parcels%get_aspect_ratio(evals)

                parcel_stats(IDX_AVG_LAM) = parcel_stats(IDX_AVG_LAM) + lam
                parcel_stats(IDX_STD_LAM) = parcel_stats(IDX_STD_LAM) + lam ** 2

                parcel_stats(IDX_SUM_VOL) = parcel_stats(IDX_SUM_VOL) + vol
                parcel_stats(IDX_STD_VOL) = parcel_stats(IDX_STD_VOL) + vol ** 2

                if (vol <= vmin) then
                    parcel_stats(IDX_N_SMALL) = parcel_stats(IDX_N_SMALL) + one
                endif

                bmax = max(bmax, b)

                bmin = min(bmin, b)

#ifndef NDEBUG
                !$omp critical
                if (abs(parcels%get_determinant(n) / parcels%get_abc(vol) ** 2 - one) > thres) then
                    call mpi_exit_on_error("Parcel determinant not preserved!")
                endif
                !$omp end critical
#endif
                parcel_stats(IDX_RMS_XI)   = parcel_stats(IDX_RMS_XI)   + vol * parcels%vorticity(1, n) ** 2
                parcel_stats(IDX_RMS_ETA)  = parcel_stats(IDX_RMS_ETA)  + vol * parcels%vorticity(2, n) ** 2
                parcel_stats(IDX_RMS_ZETA) = parcel_stats(IDX_RMS_ZETA) + vol * parcels%vorticity(3, n) ** 2

            enddo
            !$omp end do
            !$omp end parallel

            parcel_stats(IDX_MIN_BUOY) = bmin
            parcel_stats(IDX_MAX_BUOY) = bmax

            parcel_stats(IDX_NSPLITS)   = n_parcel_splits
            parcel_stats(IDX_NBIG_ICLO) = n_big_close
            parcel_stats(IDX_NMERGES)   = n_parcel_merges

            call mpi_blocking_reduce(parcel_stats(IDX_APE:IDX_NMERGES), MPI_SUM, world)
            call mpi_blocking_reduce(parcel_stats(IDX_MIN_BUOY), MPI_MIN, world)
            call mpi_blocking_reduce(parcel_stats(IDX_MAX_BUOY), MPI_MAX, world)

            parcel_merge_stats = n_way_parcel_mergers
            call mpi_blocking_reduce(parcel_merge_stats, MPI_SUM, world)

            parcels%total_num = nint(parcel_stats(IDX_NTOT_PAR), kind=int64)
            ntoti = one / dble(parcels%total_num)

            ! divide by domain volume to get domain-averaged quantities
            parcel_stats(IDX_KE) = f12 * parcel_stats(IDX_KE) * vdomaini
            if (ape_calculation == 'ape density') then
                parcel_stats(IDX_APE) = parcel_stats(IDX_APE) * vdomaini
            endif
            parcel_stats(IDX_ENSTROPHY) = f12 * parcel_stats(IDX_ENSTROPHY)* vdomaini

            parcel_stats(IDX_AVG_LAM) = parcel_stats(IDX_AVG_LAM) * ntoti
            parcel_stats(IDX_STD_LAM) = sqrt(abs(parcel_stats(IDX_STD_LAM) * ntoti &
                                                - parcel_stats(IDX_AVG_LAM) ** 2))

            parcel_stats(IDX_RMS_XI)   = sqrt(parcel_stats(IDX_RMS_XI)   / parcel_stats(IDX_SUM_VOL))
            parcel_stats(IDX_RMS_ETA)  = sqrt(parcel_stats(IDX_RMS_ETA)  / parcel_stats(IDX_SUM_VOL))
            parcel_stats(IDX_RMS_ZETA) = sqrt(parcel_stats(IDX_RMS_ZETA) / parcel_stats(IDX_SUM_VOL))

            parcel_stats(IDX_AVG_VOL) = parcel_stats(IDX_SUM_VOL) * ntoti
            parcel_stats(IDX_STD_VOL) = sqrt(abs(parcel_stats(IDX_STD_VOL) * ntoti &
                                                - parcel_stats(IDX_AVG_VOL) ** 2))

            call stop_timer(parcel_stats_timer)

        end subroutine calculate_parcel_diagnostics

end module parcel_diagnostics
