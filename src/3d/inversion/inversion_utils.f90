module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, upper, lower
    use mpi_layout
    use mpi_environment
    use sta3dfft, only : initialise_fft &
                       , finalise_fft   &
                       , rkx            &
                       , rky            &
                       , rkz            &
                       , fftxyp2s       &
                       , fftxys2p       &
                       , fftsine        &
                       , fftcosine, zfactors, ztrig
    use stafft, only : dst
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, protected, allocatable :: k2l2i(:, :)

    ! Note k2l2 = k^2+l^2
    double precision, protected, allocatable :: k2l2(:, :)

    double precision, protected, allocatable :: green(:, :, :)

    !De-aliasing filter:
    double precision, protected, allocatable :: filt(:, :, :)

    double precision, protected, allocatable :: gamtop(:), gambot(:)

    double precision, protected, allocatable :: thetam(:, :, :)    ! theta_{-}
    double precision, protected, allocatable :: thetap(:, :, :)    ! theta_{+}
    double precision, protected, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
    double precision, protected, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
    double precision, protected, allocatable :: phim(:, :, :)      ! phi_{-}
    double precision, protected, allocatable :: phip(:, :, :)      ! phi_{+}
    double precision, protected, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
    double precision, protected, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz

    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq, ap

    logical :: is_initialised = .false.

    public :: init_inversion        &
            , finalise_inversion    &
            , diffz                 &
            , central_diffz         &
            , lapinv1               &
            , dz2                   &
            , filt                  &
            , hdzi                  &
            , k2l2i                 &
            , green                 &
            , thetap                &
            , thetam                &
            , dthetap               &
            , dthetam               &
            , gambot                &
            , gamtop

    public :: field_combine_semi_spectral   &
            , field_combine_physical        &
            , field_decompose_semi_spectral &
            , field_decompose_physical      &
            , central_diffz_semi_spectral

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            double precision, allocatable :: a0(:, :)
            double precision              :: kxmaxi, kymaxi, kzmaxi
            double precision              :: skx(box%lo(1):box%hi(1)), &
                                             sky(box%lo(2):box%hi(2)), &
                                             skz(0:nz)
            integer                       :: kx, ky, iz, kz
            double precision              :: z, zm(0:nz), zp(0:nz)
            double precision              :: phip00(0:nz)

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call initialise_fft(extent)

            dz = dx(3)
            dzi = dxi(3)
            dz6  = f16 * dx(3)
            dz2  = f12 * dx(3)
            dz24 = f124 * dx(3)
            dzisq = dxi(3) ** 2
            hdzi = f12 * dxi(3)

            allocate(k2l2i(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(k2l2(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            allocate(filt(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            !Squared wavenumber array (used in tridiagonal solve):
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    k2l2(ky, kx) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                k2l2(0, 0) = one
            endif

            k2l2i = one / k2l2

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                k2l2(0, 0) = zero
                k2l2i(0, 0) = zero
            endif

            !----------------------------------------------------------
            !Define Hou and Li filter (2D and 3D):
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx(box%lo(1):box%hi(1))) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky(box%lo(2):box%hi(2))) ** 36
            kzmaxi = one/maxval(rkz)
            skz = -36.d0 * (kzmaxi * rkz) ** 36

            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                  filt(0,  ky, kx) = dexp(skx(kx) + sky(ky))
                  filt(nz, ky, kx) = filt(0, ky, kx)
                  do kz = 1, nz-1
                     filt(kz, ky, kx) = filt(0, ky, kx) * dexp(skz(kz))
                  enddo
               enddo
            enddo

            !Ensure filter does not change domain mean:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                filt(:, 0, 0) = one
            endif


            allocate(green(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            allocate(phim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(phip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dphim(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dphip(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(thetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(thetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dthetam(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(dthetap(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            !---------------------------------------------------------------------
            !Define Green function
            !$omp parallel do
            do kz = 1, nz
                green(kz, :, :) = - one / (k2l2 + rkz(kz) ** 2)
            enddo
            !$omp end parallel do
            green(0, :, :) = - k2l2i

            !---------------------------------------------------------------------
            !Define zm = zmax - z, zp = z - zmin
            !$omp parallel do private(z)
            do iz = 0, nz
                z = lower(3) + dx(3) * dble(iz)
                zm(iz) = upper(3) - z
                zp(iz) = z - lower(3)
            enddo
            !$omp end parallel do

            !Hyperbolic functions used for solutions of Laplace's equation:
            do kx = box%lo(1), box%hi(1)
                do ky = max(1, box%lo(2)), box%hi(2)
                    call set_hyperbolic_functions(kx, ky, zm, zp)
                enddo
            enddo

            ! ky = 0
            if (box%lo(2) == 0) then
                do kx = max(1, box%lo(1)), box%hi(1)
                    call set_hyperbolic_functions(kx, 0, zm, zp)
                enddo
            endif

            phip00 = zero
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                !$omp parallel workshare
                ! kx = ky = 0
                phim(:, 0, 0) = zm / extent(3)
                phip(:, 0, 0) = zp / extent(3)

                dphim(:, 0, 0) = - one / extent(3)
                dphip(:, 0, 0) =   one / extent(3)

                thetam(:, 0, 0) = zero
                thetap(:, 0, 0) = zero

                dthetam(:, 0, 0) = zero
                dthetap(:, 0, 0) = zero

                phip00 = phip(:, 0, 0)
                !$omp end parallel workshare
            endif

            !---------------------------------------------------------------------
            !Define gamtop as the integral of phip(iz, 0, 0) with zero average:

            call MPI_Allreduce(MPI_IN_PLACE,            &
                               phip00(0:nz),            &
                               nz+1,                    &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_SUM,                 &
                               world%comm,              &
                               world%err)

            !$omp parallel workshare
            gamtop = f12 * extent(3) * (phip00 ** 2 - f13)
            !$omp end parallel workshare

            !$omp parallel do
            do iz = 0, nz
                gambot(iz) = gamtop(nz-iz)
            enddo
            !$omp end parallel do
            !Here gambot is the complement of gamtop.


            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:

            allocate(a0(box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(etdv(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
            allocate(htdv(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

            a0 = -two * dzisq - k2l2
            ap = dzisq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for the vertical vorticity component:
            htdv(0, :, :) = one / a0
            etdv(0, :, :) = -two * ap * htdv(0, :, :)
            !$omp parallel shared(a0, ap, etdv, htdv, nz) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                htdv(iz, :, :) = one / (a0(:, :) + ap * etdv(iz-1, :, :))
                etdv(iz, :, :) = -ap * htdv(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                etdv(nz-1, 0, 0) = zero
            endif

            htdv(nz, :, :) = one / (a0 + two * ap * etdv(nz-1, :, :))

            ! Remove horizontally-averaged part (done separately):
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                htdv(:, 0, 0) = zero
                etdv(:, 0, 0) = zero
            endif

            deallocate(a0)

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine finalise_inversion
            deallocate(green)
            deallocate(gamtop)
            deallocate(gambot)
            deallocate(phim)
            deallocate(phip)
            deallocate(dphim)
            deallocate(dphip)
            deallocate(thetam)
            deallocate(thetap)
            deallocate(dthetam)
            deallocate(dthetap)

            deallocate(k2l2i)
            deallocate(k2l2)
            deallocate(filt)
            deallocate(etdv)
            deallocate(htdv)

            call finalise_fft

        end subroutine finalise_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_hyperbolic_functions(kx, ky, zm, zp)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zm(0:nz), zp(0:nz)
            double precision             :: R(0:nz), Q(0:nz), k2ifac
            double precision             :: ef, em(0:nz), ep(0:nz), Lm(0:nz), Lp(0:nz)
            double precision             :: fac, div, kl

            kl = dsqrt(k2l2(ky, kx))
            fac = kl * extent(3)
            ef = dexp(- fac)
#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ef is really small.
            ef = max(ef, dsqrt(tiny(ef)))
#endif
            div = one / (one - ef**2)
            k2ifac = f12 * k2l2i(ky, kx)

            Lm = kl * zm
            Lp = kl * zp

            ep = dexp(- Lp)
            em = dexp(- Lm)

#ifndef NDEBUG
            ! To avoid "Floating-point exception - erroneous arithmetic operation"
            ! when ep and em are really small.
            ep = max(ep, dsqrt(tiny(ep)))
            em = max(em, dsqrt(tiny(em)))
#endif

            phim(:, ky, kx) = div * (ep - ef * em)
            phip(:, ky, kx) = div * (em - ef * ep)

            dphim(:, ky, kx) = - kl * div * (ep + ef * em)
            dphip(:, ky, kx) =   kl * div * (em + ef * ep)

            Q = div * (one + ef**2)
            R = div * two * ef

            thetam(:, ky, kx) = k2ifac * (R * Lm * phip(:, ky, kx) - Q * Lp * phim(:, ky, kx))
            thetap(:, ky, kx) = k2ifac * (R * Lp * phim(:, ky, kx) - Q * Lm * phip(:, ky, kx))

            dthetam(:, ky, kx) = - k2ifac * ((Q * Lp - one) * dphim(:, ky, kx) - R * Lm * dphip(:, ky, kx))
            dthetap(:, ky, kx) = - k2ifac * ((Q * Lm - one) * dphip(:, ky, kx) - R * Lp * dphim(:, ky, kx))

        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! fc  - complete field (physical space)
        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! cfc - copy of complete field (physical space)
        subroutine field_decompose_physical(fc, sf)
            double precision, intent(in)  :: fc(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision, intent(out) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            call fftxyp2s(fc, sf)

            call field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : complete field (semi-spectral space)
        ! out: full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        subroutine field_decompose_semi_spectral(sfc)
            double precision, intent(inout) :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: sfctop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz, kx, ky

            ! subtract harmonic part
            !$omp parallel do
            do iz = 1, nz-1
                sfc(iz, :, :) = sfc(iz, :, :) - (sfc(0, :, :) * phim(iz, :, :) + sfc(nz, :, :) * phip(iz, :, :))
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sfctop = sfc(nz, :, :)
            !$omp end parallel workshare

            ! transform interior to fully spectral
            !$omp parallel do collapse(2)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dst(1, nz, sfc(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sfc(nz, :, :) = sfctop
            !$omp end parallel workshare

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! sf  - full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! fc  - complete field (physical space)
        ! sfc - complete field (semi-spectral space)
        subroutine field_combine_physical(sf, fc)
            double precision, intent(in)  :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: fc(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision              :: sfc(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))

            sfc = sf

            call field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! in : full-spectral (1:nz-1), semi-spectral at iz = 0 and iz = nz
        ! out: complete field (semi-spectral space)
        subroutine field_combine_semi_spectral(sf)
            double precision, intent(inout) :: sf(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: sftop(box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz, kx, ky

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
            !$omp parallel workshare
            sftop = sf(nz, :, :)
            !$omp end parallel workshare

            !$omp parallel do collapse(2)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    sf(nz, ky, kx) = zero
                    call dst(1, nz, sf(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel workshare
            sf(nz, :, :) = sftop
            !$omp end parallel workshare

            ! add harmonic part to sfc:
            !$omp parallel do
            do iz = 1, nz-1
                sf(iz, :, :) = sf(iz, :, :) + sf(0, :, :) * phim(iz, :, :) + sf(nz, :, :) * phip(iz, :, :)
            enddo
            !$omp end parallel do

        end subroutine field_combine_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f using 2nd-order differencing.
        !Here fs = f, ds = df/dz. In semi-spectral space.
        subroutine central_diffz_semi_spectral(fs, ds)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: iz

            ! Linear extrapolation at the boundaries:
            ! iz = 0:  (fs(1) - fs(0)) / dz
            ! iz = nz: (fs(nz) - fs(nz-1)) / dz
            !$omp parallel workshare
            ds(0,  :, :) = dzi * (fs(1,    :, :) - fs(0,    :, :))
            ds(nz, :, :) = dzi * (fs(nz,   :, :) - fs(nz-1, :, :))
            !$omp end parallel workshare

            ! central differencing for interior cells
            !$omp parallel do private(iz) default(shared)
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end parallel do

        end subroutine central_diffz_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f using 2nd-order differencing.
        !Here df = df/dz. In physical space.
        subroutine central_diffz(f, df)
            double precision, intent(in)  :: f(-1:nz+1,  box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision, intent(out) :: df(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            integer                       :: iz

            ! Linear extrapolation at the boundaries:
            ! iz = 0:  (f(1) - f(0)) / dz
            ! iz = nz: (f(nz) - f(nz-1)) / dz
            !$omp parallel workshare
            df(0,  :, :) = dzi * (f(1,    :, :) - f(0,    :, :))
            df(nz, :, :) = dzi * (f(nz,   :, :) - f(nz-1, :, :))
            !$omp end parallel workshare

            ! central differencing for interior cells
            !$omp parallel do private(iz) default(shared)
            do iz = 1, nz-1
                df(iz, :, :) = (f(iz+1, :, :) - f(iz-1, :, :)) * hdzi
            enddo
            !$omp end parallel do

            ! Linear extrapolate for halo grid points:
            !$omp parallel workshare
            df(  -1, :, :) = two * df( 0, :, :) - df(   1, :, :)
            df(nz+1, :, :) = two * df(nz, :, :) - df(nz-1, :, :)
            !$omp end parallel workshare

        end subroutine central_diffz

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f in mixed-spectral space
        !Here fs = f, ds = df/dz. Both fields are in mixed-spectral space.
        ! fs - mixed-spectral space
        ! ds - derivative linear part
        ! as - derivative sine part
        subroutine diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision              :: as(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                       :: kz, iz

            !Calculate the derivative of the linear part (ds) in semi-spectral space:
            !$omp parallel do private(iz)  default(shared)
            do iz = 0, nz
                ds(iz, :, :) = fs(0, :, :) * dphim(iz, :, :) + fs(nz, :, :) * dphip(iz, :, :)
            enddo
            !$omp end parallel do

            ! Calculate d/dz of this sine series:
            !$omp parallel workshare
            as(0, :, :) = zero
            !$omp end parallel workshare
            !$omp parallel do private(kz)  default(shared)
            do kz = 1, nz-1
                as(kz, :, :) = rkz(kz) * fs(kz, :, :)
            enddo
            !$omp end parallel do
            !$omp parallel workshare
            as(nz, :, :) = zero
            !$omp end parallel workshare

            !FFT these quantities back to semi-spectral space:
            call fftcosine(as)

            ! Combine vertical derivative given the sine (as) and linear (ds) parts:
            !omp parallel workshare
            ds = ds + as
            !omp end parallel workshare

            call field_decompose_semi_spectral(ds)

        end subroutine diffz

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here dfs/dz = 0 on the z boundaries.
        !Uses 2nd-order differencing
        !*** Overwrites fs ***
        subroutine lapinv1(fs)
            double precision, intent(inout) :: fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            double precision                :: rs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1))
            integer                         :: iz

            rs = fs
            fs(0, :, :) = rs(0, :, :) * htdv(0, :, :)

            !$omp parallel shared(rs, fs, ap, htdv, nz) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                fs(iz, :, :) = (rs(iz, :, :) - ap * fs(iz-1, :, :)) * htdv(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = (rs(nz, :, :) - two * ap * fs(nz-1, :, :)) * htdv(nz, :, :)

            !$omp parallel shared(fs, etdv, nz) private(iz) default(none)
            !$omp do
            do iz = nz-1, 0, -1
                fs(iz, :, :) = etdv(iz, :, :) * fs(iz+1, :, :) + fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Zero horizontal wavenumber in x & y treated separately:
            if ((box%lo(1) == 0) .and. (box%lo(2) == 0)) then
                fs(:, 0, 0) = zero
            endif
        end subroutine

end module inversion_utils
