module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, upper, lower
    use mpi_layout
    use mpi_communicator
    use sta3dfft, only : initialise_fft &
                       , finalise_fft   &
                       , rkx            &
                       , rky            &
                       , rkz            &
                       , fftxyp2s       &
                       , fftxys2p
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, allocatable :: k2l2i(:, :)

    ! Note k2l2 = k^2+l^2
    double precision, allocatable :: k2l2(:, :)

    integer, parameter :: nsubs_tri = 8 !number of blocks for openmp
    integer :: nxsub

    double precision, allocatable :: green(:, :, :)

    !De-aliasing filter:
    double precision, allocatable :: filt(:, :, :)

    double precision, allocatable :: gamtop(:), gambot(:)

    double precision, allocatable :: thetam(:, :, :)    ! theta_{-}
    double precision, allocatable :: thetap(:, :, :)    ! theta_{+}
    double precision, allocatable :: dthetam(:, :, :)   ! dtheta_{-}/dz
    double precision, allocatable :: dthetap(:, :, :)   ! dtheta_{+}/dz
    double precision, allocatable :: phim(:, :, :)      ! phi_{-}
    double precision, allocatable :: phip(:, :, :)      ! phi_{+}
    double precision, allocatable :: dphim(:, :, :)     ! dphi_{-}/dz
    double precision, allocatable :: dphip(:, :, :)     ! dphi_{+}/dz

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
            , field_decompose_physical

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            double precision, allocatable :: a0(:, :)
            double precision              :: kxmaxi, kymaxi, kzmaxi
            double precision              :: skx(0:nx-1), sky(0:ny-1), skz(0:nz)
            integer                       :: isub, ib_sub, ie_sub
            integer                       :: kx, ky, iz, kz
            double precision              :: z, zm(0:nz), zp(0:nz)

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call initialise_fft(extent)

            allocate(green(0:nz, 0:nx-1, 0:ny-1))
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            allocate(phim(0:nz, 0:nx-1, 0:ny-1))
            allocate(phip(0:nz, 0:nx-1, 0:ny-1))
            allocate(dphim(0:nz, 0:nx-1, 0:ny-1))
            allocate(dphip(0:nz, 0:nx-1, 0:ny-1))
            allocate(thetam(0:nz, 0:nx-1, 0:ny-1))
            allocate(thetap(0:nz, 0:nx-1, 0:ny-1))
            allocate(dthetam(0:nz, 0:nx-1, 0:ny-1))
            allocate(dthetap(0:nz, 0:nx-1, 0:ny-1))

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
            do ky = 1, ny-1
                do kx = 0, nx-1
                    call set_hyperbolic_functions(kx, ky, zm, zp)
                enddo
            enddo

            ! ky = 0
            do kx = 1, nx-1
                call set_hyperbolic_functions(kx, 0, zm, zp)
            enddo

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
            !$omp end parallel workshare

            !---------------------------------------------------------------------
            !Define gamtop as the integral of phip(iz, 0, 0) with zero average:
            !$omp parallel workshare
            gamtop = f12 * extent(3) * (phip(:, 0, 0) ** 2 - f13)
            !$omp end parallel workshare

            !$omp parallel do
            do iz = 0, nz
                gambot(iz) = gamtop(nz-iz)
            enddo
            !$omp end parallel do
            !Here gambot is the complement of gamtop.



            dz = dx(3)
            dzi = dxi(3)
            dz6  = f16 * dx(3)
            dz2  = f12 * dx(3)
            dz24 = f124 * dx(3)
            dzisq = dxi(3) ** 2
            hdzi = f12 * dxi(3)

            allocate(k2l2i(0:nx-1, 0:ny-1))
            allocate(k2l2(0:nx-1, 0:ny-1))

            allocate(filt(0:nz, 0:nx-1, 0:ny-1))

            allocate(a0(nx, ny))
            allocate(etdv(0:nz, nx, ny))
            allocate(htdv(0:nz, nx, ny))

            nxsub = nx / nsubs_tri

            !Squared wavenumber array (used in tridiagonal solve):
            do ky = 0, ny-1
                do kx = 0, nx-1
                    k2l2(kx, ky) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            k2l2(0, 0) = one
            k2l2i = one / k2l2
            k2l2(0, 0) = zero
            k2l2i(0, 0) = zero

            !----------------------------------------------------------
            !Define Hou and Li filter (2D and 3D):
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky) ** 36
            kzmaxi = one/maxval(rkz)
            skz = -36.d0 * (kzmaxi * rkz) ** 36

            do ky = 0, ny-1
               do kx = 0, nx-1
                     !filt(:, kx, ky) = dexp(skx(kx) + sky(ky))
                  filt(0,  kx, ky) = dexp(skx(kx) + sky(ky))
                  filt(nz, kx, ky) = filt(0, kx, ky)
                  do kz = 1, nz-1
                     filt(kz, kx, ky) = filt(0, kx, ky) * dexp(skz(kz))
                  enddo
               enddo
            enddo

            !Ensure filter does not change domain mean:
            filt(:, 0, 0) = one

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0 = -two * dzisq - k2l2
            ap = dzisq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for the vertical vorticity component:
            htdv(0, :, :) = one / a0
            etdv(0, :, :) = -two * ap * htdv(0, :, :)
            !$omp parallel shared(a0, ap, etdv, htdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    htdv(iz, ib_sub:ie_sub, :) = one / (a0(ib_sub:ie_sub, :) &
                                               + ap * etdv(iz-1, ib_sub:ie_sub, :))
                    etdv(iz, ib_sub:ie_sub, :) = -ap * htdv(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            etdv(nz-1, 1, 1) = zero

            htdv(nz, :, :) = one / (a0 + two * ap * etdv(nz-1, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdv(:, 1, 1) = zero
            etdv(:, 1, 1) = zero

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

            kl = dsqrt(k2l2(kx, ky))
            fac = kl * extent(3)
            ef = dexp(- fac)
            div = one / (one - ef**2)
            k2ifac = f12 * k2l2i(kx, ky)

            Lm = kl * zm
            Lp = kl * zp

            ep = dexp(- Lp)
            em = dexp(- Lm)

            phim(:, kx, ky) = div * (ep - ef * em)
            phip(:, kx, ky) = div * (em - ef * ep)

            dphim(:, kx, ky) = - kl * div * (ep + ef * em)
            dphip(:, kx, ky) =   kl * div * (em + ef * ep)

            Q = div * (one + ef**2)
            R = div * two * ef

            thetam(:, kx, ky) = k2ifac * (R * Lm * phip(:, kx, ky) - Q * Lp * phim(:, kx, ky))
            thetap(:, kx, ky) = k2ifac * (R * Lp * phim(:, kx, ky) - Q * Lm * phip(:, kx, ky))

            dthetam(:, kx, ky) = - k2ifac * ((Q * Lp - one) * dphim(:, kx, ky) - R * Lm * dphip(:, kx, ky))
            dthetap(:, kx, ky) = - k2ifac * ((Q * Lm - one) * dphip(:, kx, ky) - R * Lp * dphim(:, kx, ky))

        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_decompose_physical(fc, sf)
            double precision, intent(in)  :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision, intent(out) :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision              :: cfc(0:nz, 0:ny-1, 0:nx-1)   ! copy of complete field (physical space)

            cfc = fc
            call fftxyp2s(cfc, sf)

            call field_decompose_semi_spectral(sf)

        end subroutine field_decompose_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_decompose_semi_spectral(sfc)
            double precision, intent(inout) :: sfc(0:nz, 0:nx-1, 0:ny-1) ! in : complete field (semi-spectral space)
                                                                         ! out: full-spectral (1:nz-1),
                                                                         !      semi-spectral at iz = 0 and iz = nz
            double precision                :: sfctop(0:nx-1, 0:ny-1)
            integer                         :: iz

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
            !call fftsine(sfc) !FIXME enable

            !$omp parallel workshare
            sfc(nz, :, :) = sfctop
            !$omp end parallel workshare

        end subroutine field_decompose_semi_spectral

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_combine_physical(sf, fc)
            double precision, intent(in)  :: sf(0:nz, 0:nx-1, 0:ny-1)    ! full-spectral (1:nz-1),
                                                                         ! semi-spectral at iz = 0 and iz = nz
            double precision, intent(out) :: fc(0:nz, 0:ny-1, 0:nx-1)    ! complete field (physical space)
            double precision              :: sfc(0:nz, 0:nx-1, 0:ny-1)   ! complete field (semi-spectral space)

            sfc = sf

            call field_combine_semi_spectral(sfc)

            ! transform to physical space as fc:
            call fftxys2p(sfc, fc)

        end subroutine field_combine_physical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine field_combine_semi_spectral(sf)
            double precision, intent(inout) :: sf(0:nz, 0:nx-1, 0:ny-1) ! in: full-spectral (1:nz-1),
                                                                        !     semi-spectral at iz = 0 and iz = nz
                                                                        ! out: complete field (semi-spectral space)
            double precision                :: sftop(0:nx-1, 0:ny-1)
            integer                         :: iz

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
            !$omp parallel workshare
            sftop = sf(nz, :, :)
            !$omp end parallel workshare

            ! call fftsine(sf) !FIXME enable

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
        !Here fs = f, ds = df/dz. In semi-spectral space or physical space.
        subroutine central_diffz(fs, ds)
            double precision, intent(in)  :: fs(0:, 0:, 0:)
            double precision, intent(out) :: ds(0:, 0:, 0:)
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

        end subroutine central_diffz

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f in mixed-spectral space
        !Here fs = f, ds = df/dz. Both fields are in mixed-spectral space.
        subroutine diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, 0:nx-1, 0:ny-1) ! f in mixed-spectral space
            double precision, intent(out) :: ds(0:nz, 0:nx-1, 0:ny-1) ! derivative linear part
            double precision              :: as(0:nz, 0:nx-1, 0:ny-1) ! derivative sine-part
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
!             call fftcosine(as) !FIXME enable

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
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            double precision                :: rs(0:nz, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            rs = fs
            fs(0, :, :) = rs(0, :, :) * htdv(0, :, :)

            !$omp parallel shared(rs, fs, ap, htdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    fs(iz, ib_sub:ie_sub, :) = (rs(iz, ib_sub:ie_sub, :) &
                                             - ap * fs(iz-1, ib_sub:ie_sub, :)) &
                                             * htdv(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = (rs(nz, :, :) - two * ap * fs(nz-1, :, :)) * htdv(nz, :, :)

            !$omp parallel shared(fs, etdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-1, 0, -1
                    fs(iz, ib_sub:ie_sub, :) = etdv(iz, ib_sub:ie_sub, :) &
                                             * fs(iz+1, ib_sub:ie_sub, :) + fs(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

             !Zero horizontal wavenumber in x & y treated separately:
             fs(:, 1, 1) = zero
        end subroutine

end module inversion_utils
