! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse
    use h5_utils
    use h5_writer
    use omp_lib
    implicit none

    private

    ! peref : potential energy reference
    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: peref, pe, ke

#ifdef ENABLE_DIAGNOSE
    ! avg_lam : mean aspect ratio over all parcels
    ! avg_vol : mean volume over all parcels
    ! std_lam : standard deviation of aspect ratio
    ! std_vol : standard deviation of volume
    double precision :: avg_lam, avg_vol
    double precision :: std_lam, std_vol

    ! buoyancy weighted first and second moments
    double precision :: xb_bar, x2b_bar
    double precision :: zb_bar, z2b_bar
    double precision :: xzb_bar

    ! vorticity weighted first and second moments
    double precision :: xv_bar, x2v_bar
    double precision :: zv_bar, z2v_bar
    double precision :: xzv_bar
#endif

    public :: init_parcel_diagnostics, &
              write_h5_parcel_diagnostics

    contains

        ! compute the reference potential energy
        subroutine init_parcel_diagnostics
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = one / extent(1)
            zmean = lower(2) + f12 * gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * parcels%volume(ii(n))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo
        end subroutine init_parcel_diagnostics


        subroutine calculate_diagnostics
            integer          :: n
            double precision :: b, z, vel(2), vol, zmin
#ifdef ENABLE_DIAGNOSE
            double precision :: eval, lam, B22, lsum, l2sum, vsum, v2sum
            double precision :: xbv, x2bv, zbv, z2bv, xzbv
            double precision :: xvv, x2vv, zvv, z2vv, xzvv
            double precision :: bvsum, vvsum, bv, vv
#endif
            ! reset
            ke = zero
            pe = zero

            zmin = lower(2)

#ifdef ENABLE_DIAGNOSE
            avg_lam = zero
            avg_vol = zero
            std_lam = zero
            std_vol = zero

            xb_bar = zero
            zb_bar = zero
            x2b_bar = zero
            z2b_bar = zero

            xv_bar = zero
            zv_bar = zero
            x2v_bar = zero
            z2v_bar = zero

            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z, eval, lam, B22, bv, vv) &
            !$omp& reduction(+: ke, pe, lsum, l2sum, vsum, v2sum, bvsum, xbv, &
            !$omp&              zbv, x2bv, z2bv, xzbv, xvv, zvv, x2vv, z2vv, xzvv)
#else
            !$omp parallel default(shared)
            !$omp do private(n, vel, vol, b, z) reduction(+: ke, pe)
#endif
            do n = 1, n_parcels

                vel = parcels%velocity(n, :)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(n, 2) - zmin

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                ! potential energy
                pe = pe - b * z * vol

#ifdef ENABLE_DIAGNOSE
                B22 = get_B22(parcels%B(n, 1), parcels%B(n, 2), vol)
                eval = get_eigenvalue(parcels%B(n, 1), parcels%B(n, 2), B22)
                lam = get_aspect_ratio(eval, vol)

                lsum = lsum + lam
                l2sum = l2sum + lam ** 2

                vsum = vsum + vol
                v2sum = v2sum + vol ** 2

                ! we only use the upper half in horizontal direction
                if (parcels%position(n, 1) >= 0) then
                    bv = parcels%buoyancy(n) * vol
                    bvsum = bvsum + bv
                    xbv = xbv + bv * parcels%position(n, 1)
                    zbv = zbv + bv * parcels%position(n, 2)

                    x2bv = x2bv + bv * parcels%position(n, 1) ** 2
                    z2bv = z2bv + bv * parcels%position(n, 2) ** 2
                    xzbv = xzbv + bv * parcels%position(n, 2) * parcels%position(n, 2)


                    vv = parcels%vorticity(n) * vol
                    vvsum = vvsum + vv
                    xvv = xvv + vv * parcels%position(n, 1)
                    zvv = zvv + vv * parcels%position(n, 2)

                    x2vv = x2vv + vv * parcels%position(n, 1) ** 2
                    z2vv = z2vv + vv * parcels%position(n, 2) ** 2
                    xzvv = xzvv + vv * parcels%position(n, 2) * parcels%position(n, 2)
                endif
#endif
            enddo
            !$omp end do
            !$omp end parallel

            ke = f12 * ke
            pe = pe - peref

#ifdef ENABLE_DIAGNOSE
            avg_lam = lsum / dble(n_parcels)
            std_lam = dsqrt(l2sum / dble(n_parcels) - avg_vol ** 2)

            avg_vol = vsum / dble(n_parcels)
            std_vol = dsqrt(v2sum / dble(n_parcels) - avg_vol ** 2)

            ! we do not need to divide by the number of of involved
            ! parcels since whe divide by the sums "bvsum" or "vvsum"
            ! what should also be averages (i.e. divided by the number of
            ! involved parcels)
            bvsum = one / bvsum
            vvsum = one / vvsum

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
#endif
        end subroutine calculate_diagnostics


        subroutine write_h5_parcel_diagnostics(h5file_id)
            integer(hid_t), intent(in)    :: h5file_id

            call calculate_diagnostics

            !
            ! write diagnostics
            !
            call write_h5_double_scalar_attrib(h5file_id, "potential energy", pe)
            call write_h5_double_scalar_attrib(h5file_id, "kinetic energy", ke)
            call write_h5_double_scalar_attrib(h5file_id, "total energy", ke + pe)
            call write_h5_int_scalar_attrib(h5file_id, "num parcel", n_parcels)

#ifdef ENABLE_DIAGNOSE
            call write_h5_double_scalar_attrib(h5file_id, "avg aspect ratio", avg_lam)
            call write_h5_double_scalar_attrib(h5file_id, "std aspect ratio", std_lam)
            call write_h5_double_scalar_attrib(h5file_id, "avg volume", avg_vol)
            call write_h5_double_scalar_attrib(h5file_id, "std volume", std_vol)

            call write_h5_double_scalar_attrib(h5file_id, "xb_bar", xb_bar)
            call write_h5_double_scalar_attrib(h5file_id, "x2b_bar", x2b_bar)
            call write_h5_double_scalar_attrib(h5file_id, "zb_bar", zb_bar)
            call write_h5_double_scalar_attrib(h5file_id, "z2b_bar", z2b_bar)
            call write_h5_double_scalar_attrib(h5file_id, "xzb_bar", xzb_bar)

            call write_h5_double_scalar_attrib(h5file_id, "xv_bar", xv_bar)
            call write_h5_double_scalar_attrib(h5file_id, "x2v_bar", x2v_bar)
            call write_h5_double_scalar_attrib(h5file_id, "zv_bar", zv_bar)
            call write_h5_double_scalar_attrib(h5file_id, "z2v_bar", z2v_bar)
            call write_h5_double_scalar_attrib(h5file_id, "xzv_bar", xzv_bar)
#endif
        end subroutine write_h5_parcel_diagnostics
end module parcel_diagnostics
