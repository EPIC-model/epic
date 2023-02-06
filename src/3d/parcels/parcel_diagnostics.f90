! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12, thousand
    use merge_sort
    use parameters, only : extent, lower, vcell, vmin, nx, nz
    use parcel_container, only : parcels, n_parcels, n_total_parcels
    use parcel_ellipsoid
    use parcel_split_mod, only : n_parcel_splits
    use parcel_merge, only : n_parcel_merges
    use omp_lib
    use physics, only : peref
    use mpi_timer, only : start_timer, stop_timer
    use mpi_communicator
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    integer :: parcel_stats_timer

#ifndef NDEBUG
    double precision, parameter :: thres = thousand * epsilon(zero)
#endif

    integer, parameter :: IDX_PE        =  1, & ! potential energy
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
                          IDX_NMERGES   = 15, & ! number of parcel merges since last write
                          IDX_MIN_BUOY  = 16, & ! minimum parcel buoyancy
                          IDX_MAX_BUOY  = 17    ! maximum parcel buoyancy

    double precision :: parcel_stats(IDX_MAX_BUOY)

    ! buoyancy extrema
    double precision :: bmin, bmax

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

            gam = one / (extent(1) * extent(2))

            ! initially all parcels have the same volume,
            ! hence, this is valid on each MPI rank
            zmean = gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * (parcels%volume(ii(n-1)) + parcels%volume(ii(n)))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo

            call mpi_blocking_reduce(peref, MPI_SUM)

            call stop_timer(parcel_stats_timer)
        end subroutine calculate_peref


        ! Calculate all parcel related diagnostics
        subroutine calculate_parcel_diagnostics(velocity)
            double precision, intent(in) :: velocity(:, :)
            integer          :: n
            double precision :: b, z, vel(3), vol, zmin, vor(3)
            double precision :: ntoti, bmin, bmax
            double precision :: evals(3), lam

            call start_timer(parcel_stats_timer)

            ! reset
            parcel_stats = zero
            bmin = zero
            bmax = zero

            zmin = lower(3)

            parcel_stats(IDX_NTOT_PAR) = dble(n_parcels)

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z, evals, lam, vor) &
            !$omp& reduction(+: parcel_stats)
            !$omp& reduction(max: bmax) reduction(min: bmin)
            do n = 1, n_parcels

                vel = velocity(:, n)
                vor = parcels%vorticity(:, n)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(3, n) - zmin

                ! kinetic energy
                parcel_stats(IDX_KE) = parcel_stats(IDX_KE) + (vel(1) ** 2 + vel(2) ** 2 + vel(3) ** 2) * vol

                ! potential energy
                parcel_stats(IDX_PE) = parcel_stats(IDX_PE) - b * z * vol

                ! enstrophy
                parcel_stats(IDX_ENSTROPHY) = parcel_stats(IDX_ENSTROPHY) &
                                            + (vor(1) ** 2 + vor(2) ** 2 + vor(3) ** 2) * vol

                evals = get_eigenvalues(parcels%B(:, n), parcels%volume(n))
                lam = get_aspect_ratio(evals)

                parcel_stats(IDX_AVG_LAM) = parcel_stats(IDX_AVG_LAM) + lam
                parcel_stats(IDX_STD_LAM) = parcel_stats(IDX_STD_LAM) + lam ** 2

                parcel_stats(IDX_SUM_VOL) = parcel_stats(IDX_SUM_VOL) + vol
                parcel_stats(IDX_STD_VOL) = parcel_stats(IDX_STD_VOL) + vol ** 2

                if (vol <= vmin) then
                    parcel_stats(IDX_N_SMALL) = parcel_stats(IDX_N_SMALL) + one
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
                parcel_stats(IDX_RMS_XI)   = parcel_stats(IDX_RMS_XI)   + vol * parcels%vorticity(1, n) ** 2
                parcel_stats(IDX_RMS_ETA)  = parcel_stats(IDX_RMS_ETA)  + vol * parcels%vorticity(2, n) ** 2
                parcel_stats(IDX_RMS_ZETA) = parcel_stats(IDX_RMS_ZETA) + vol * parcels%vorticity(3, n) ** 2

            enddo
            !$omp end do
            !$omp end parallel

            parcel_stats(IDX_MIN_BUOY) = bmin
            parcel_stats(IDX_MAX_BUOY) = bmax

            parcel_stats(IDX_NSPLITS) = n_parcel_splits
            parcel_stats(IDX_NMERGES) = n_parcel_merges

            call mpi_blocking_reduce(parcel_stats(IDX_PE:IDX_NMERGES), MPI_SUM)
            call mpi_blocking_reduce(parcel_stats(IDX_MIN_BUOY), MPI_MIN)
            call mpi_blocking_reduce(parcel_stats(IDX_MAX_BUOY), MPI_MAX)

            n_total_parcels = int(parcel_stats(IDX_NTOT_PAR))
            ntoti = one / parcel_stats(IDX_NTOT_PAR)

            parcel_stats(IDX_KE) = f12 * parcel_stats(IDX_KE)
            parcel_stats(IDX_PE) = parcel_stats(IDX_PE) - peref
            parcel_stats(IDX_ENSTROPHY) = f12 * parcel_stats(IDX_ENSTROPHY)

            parcel_stats(IDX_AVG_LAM) = parcel_stats(IDX_AVG_LAM) * ntoti
            parcel_stats(IDX_STD_LAM) = dsqrt(abs(parcel_stats(IDX_STD_LAM) * ntoti &
                                                - parcel_stats(IDX_AVG_LAM) ** 2))

            parcel_stats(IDX_RMS_XI)   = dsqrt(parcel_stats(IDX_RMS_XI)   / parcel_stats(IDX_SUM_VOL))
            parcel_stats(IDX_RMS_ETA)  = dsqrt(parcel_stats(IDX_RMS_ETA)  / parcel_stats(IDX_SUM_VOL))
            parcel_stats(IDX_RMS_ZETA) = dsqrt(parcel_stats(IDX_RMS_ZETA) / parcel_stats(IDX_SUM_VOL))

            parcel_stats(IDX_AVG_VOL) = parcel_stats(IDX_SUM_VOL) * ntoti
            parcel_stats(IDX_STD_VOL) = dsqrt(abs(parcel_stats(IDX_STD_VOL) * ntoti &
                                                - parcel_stats(IDX_AVG_VOL) ** 2))

            call stop_timer(parcel_stats_timer)

        end subroutine calculate_parcel_diagnostics

end module parcel_diagnostics
