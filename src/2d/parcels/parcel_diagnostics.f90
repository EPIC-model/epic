! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower, vcell, vmin, nx, nz
    use options, only : verbose, write_h5_options
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse
    use h5_utils
    use h5_writer
    use omp_lib
    use timer, only : start_timer, stop_timer
    implicit none


    private

    ! h5 file handle
    integer(hid_t)     :: h5file_id
    character(len=512) :: h5fname

    ! peref : potential energy reference
    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: peref, pe, ke

    integer :: n_small

    ! avg_lam : mean aspect ratio over all parcels
    ! avg_vol : mean volume over all parcels
    ! std_lam : standard deviation of aspect ratio
    ! std_vol : standard deviation of volume
    double precision :: avg_lam, avg_vol
    double precision :: std_lam, std_vol

    ! rms vorticity
    double precision :: rms_zeta

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

    integer :: hdf5_parcel_stat_timer

    public :: create_h5_parcel_stat_file,   &
              init_parcel_diagnostics,      &
              calc_parcel_diagnostics,      &
              write_h5_parcel_stats_step,   &
              hdf5_parcel_stat_timer

    contains

        ! Create the parcel diagnostic file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_h5_parcel_stat_file(basename, overwrite)
            character(*), intent(in) :: basename
            logical,      intent(in) :: overwrite

            h5fname =  basename // '_parcel_stats.hdf5'

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_scalar_attrib(h5file_id, 'output_type', 'parcel diagnostics')

            call write_h5_timestamp(h5file_id)
            call write_h5_options(h5file_id)
            call write_h5_box(h5file_id, lower, extent, (/nx, nz/))

            call close_h5_file(h5file_id)

        end subroutine create_h5_parcel_stat_file

        ! Compute the reference potential energy
        subroutine init_parcel_diagnostics
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = one / extent(1)
            zmean = f12 * gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * parcels%volume(ii(n))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo
        end subroutine init_parcel_diagnostics


        ! Calculate all parcel related diagnostics
        subroutine calc_parcel_diagnostics(velocity)
            double precision :: velocity(:, :)
            integer          :: n
            double precision :: b, z, vel(2), vol, zmin
            double precision :: eval, lam, B22, lsum, l2sum, vsum, v2sum

            ! reset
            ke = zero
            pe = zero

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
            !$omp& reduction(+: ke, pe, lsum, l2sum, vsum, v2sum, n_small, rms_zeta)
            do n = 1, n_parcels

                vel = velocity(n, :)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(n, 2) - zmin

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                ! potential energy
                pe = pe - b * z * vol

                B22 = get_B22(parcels%B(n, 1), parcels%B(n, 2), vol)
                eval = get_eigenvalue(parcels%B(n, 1), parcels%B(n, 2), B22)
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

            ke = f12 * ke
            pe = pe - peref

            avg_lam = lsum / dble(n_parcels)
            std_lam = dsqrt(abs(l2sum / dble(n_parcels) - avg_lam ** 2))

            rms_zeta = dsqrt(rms_zeta / vsum)

            avg_vol = vsum / dble(n_parcels)
            std_vol = dsqrt(abs(v2sum / dble(n_parcels) - avg_vol ** 2))


#ifdef ENABLE_DIAGNOSE
            call straka_diagnostics
#endif
        end subroutine calc_parcel_diagnostics


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
                if (parcels%position(n, 1) >= 0) then
                    bv = parcels%buoyancy(n) * parcels%volume(n)
                    bvsum = bvsum + bv
                    xbv = xbv + bv * parcels%position(n, 1)
                    zbv = zbv + bv * parcels%position(n, 2)

                    x2bv = x2bv + bv * parcels%position(n, 1) ** 2
                    z2bv = z2bv + bv * parcels%position(n, 2) ** 2
                    xzbv = xzbv + bv * parcels%position(n, 1) * parcels%position(n, 2)


                    vv = parcels%vorticity(n) * parcels%volume(n)
                    vvsum = vvsum + vv
                    xvv = xvv + vv * parcels%position(n, 1)
                    zvv = zvv + vv * parcels%position(n, 2)

                    x2vv = x2vv + vv * parcels%position(n, 1) ** 2
                    z2vv = z2vv + vv * parcels%position(n, 2) ** 2
                    xzvv = xzvv + vv * parcels%position(n, 1) * parcels%position(n, 2)
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

        ! Write a step in the parcel diagnostic file.
        ! @param[inout] nw counts the number of writes
        ! @param[in] t is the time
        ! @param[in] dt is the time step
        subroutine write_h5_parcel_stats_step(nw, t, dt)
            integer,          intent(inout) :: nw
            double precision, intent(in)    :: t
            double precision, intent(in)    :: dt
            integer(hid_t)                  :: group
            character(:), allocatable       :: name
            logical                         :: created

            call start_timer(hdf5_parcel_stat_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a19)", "write parcel diagnostics to h5"
            endif
#endif

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5file_id)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif


            call write_h5_scalar_attrib(group, "t", t)

            call write_h5_scalar_attrib(group, "dt", dt)

            !
            ! write diagnostics
            !
            call write_h5_scalar_attrib(group, "potential energy", pe)
            call write_h5_scalar_attrib(group, "kinetic energy", ke)
            call write_h5_scalar_attrib(group, "total energy", ke + pe)
            call write_h5_scalar_attrib(group, "num parcel", n_parcels)
            call write_h5_scalar_attrib(group, "num small parcels", n_small)


            call write_h5_scalar_attrib(group, "avg aspect ratio", avg_lam)
            call write_h5_scalar_attrib(group, "std aspect ratio", std_lam)
            call write_h5_scalar_attrib(group, "avg volume", avg_vol)
            call write_h5_scalar_attrib(group, "std volume", std_vol)

            call write_h5_scalar_attrib(group, "rms vorticity", rms_zeta)

#ifdef ENABLE_DIAGNOSE
            call write_h5_scalar_attrib(group, "xb_bar", xb_bar)
            call write_h5_scalar_attrib(group, "x2b_bar", x2b_bar)
            call write_h5_scalar_attrib(group, "zb_bar", zb_bar)
            call write_h5_scalar_attrib(group, "z2b_bar", z2b_bar)
            call write_h5_scalar_attrib(group, "xzb_bar", xzb_bar)

            call write_h5_scalar_attrib(group, "xv_bar", xv_bar)
            call write_h5_scalar_attrib(group, "x2v_bar", x2v_bar)
            call write_h5_scalar_attrib(group, "zv_bar", zv_bar)
            call write_h5_scalar_attrib(group, "z2v_bar", z2v_bar)
            call write_h5_scalar_attrib(group, "xzv_bar", xzv_bar)
#endif
            call close_h5_group(group)

            ! increment counter
            nw = nw + 1

            ! update number of iterations to h5 file
            call write_h5_num_steps(h5file_id, nw)


            call close_h5_file(h5file_id)

            call stop_timer(hdf5_parcel_stat_timer)

        end subroutine write_h5_parcel_stats_step
end module parcel_diagnostics
