module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent, upper, lower
    use stafft
    use sta2dfft
    use deriv1d, only : init_deriv
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    ! Wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:), rkzi(:)

    ! Note k2l2i = 1/(k^2+l^2) (except k = l = 0, then k2l2i(0, 0) = 0)
    double precision, allocatable :: k2l2i(:, :)

    ! Note k2l2 = k^2+l^2
    double precision, allocatable :: k2l2(:, :)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer :: xfactors(5), yfactors(5), zfactors(5)
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

    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq, ap
    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_initialised = .false.
    logical :: is_fft_initialised = .false.

    public :: init_inversion  &
            , init_fft        &
            , diffx           &
            , diffy           &
            , diffz           &
            , lapinv1         &
            , fftxyp2s        &
            , fftxys2p        &
            , dz2             &
            , filt            &
            , hdzi            &
            , xfactors        &
            , yfactors        &
            , zfactors        &
            , xtrig           &
            , ytrig           &
            , ztrig           &
            , rkx             &
            , rky             &
            , rkz             &
            , k2l2i           &
            , green           &
            , rkzi            &
            , thetap          &
            , thetam          &
            , dthetap         &
            , dthetam         &
            , gambot          &
            , gamtop

    public :: field_combine_semi_spectral   &
            , field_combine_physical        &
            , field_decompose_semi_spectral &
            , field_decompose_physical

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_inversion
            integer          :: kx, ky, iz, kz
            double precision :: z, zm(0:nz), zp(0:nz)

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            call init_fft

            allocate(green(0:nz, 0:nx-1, 0:ny-1))
            allocate(gamtop(0:nz))
            allocate(gambot(0:nz))

            allocate(phim(0:nz, 0:nx-1, 0:ny-1))
            allocate(phip(0:nz, 0:nx-1, 0:ny-1))
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

        end subroutine init_inversion

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! for kx > 0 and ky >= 0 or kx >= 0 and ky > 0
        subroutine set_hyperbolic_functions(kx, ky, zm, zp)
            integer,          intent(in) :: kx, ky
            double precision, intent(in) :: zm(0:nz), zp(0:nz)
            double precision             :: R(0:nz), Q(0:nz), k2ifac, dphim(0:nz), dphip(0:nz)
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

            dphim = - kl * div * (ep + ef * em)
            dphip =   kl * div * (em + ef * ep)

            Q = div * (one + ef**2)
            R = div * two * ef

            thetam(:, kx, ky) = k2ifac * (R * Lm * phip(:, kx, ky) - Q * Lp * phim(:, kx, ky))
            thetap(:, kx, ky) = k2ifac * (R * Lp * phim(:, kx, ky) - Q * Lm * phip(:, kx, ky))

            dthetam(:, kx, ky) = - k2ifac * ((Q * Lp - one) * dphim - R * Lm * dphip)
            dthetap(:, kx, ky) = - k2ifac * ((Q * Lm - one) * dphip - R * Lp * dphim)

        end subroutine set_hyperbolic_functions

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision, allocatable :: a0(:, :)
            double precision              :: rkxmax, rkymax, rksqmax
            double precision              :: kxmaxi, kymaxi, kzmaxi
            integer                       :: kx, ky, kz
            double precision              :: skx(0:nx-1), sky(0:ny-1), skz(0:nz)
            integer                       :: iz, isub, ib_sub, ie_sub

            if (is_fft_initialised) then
                return
            endif

            is_fft_initialised = .true.

            dz = dx(3)
            dzi = dxi(3)
            dz6  = f16 * dx(3)
            dz2  = f12 * dx(3)
            dz24 = f124 * dx(3)
            dzisq = dxi(3) ** 2
            hdzi = f12 * dxi(3)
            nwx = nx / 2
            nwy = ny / 2
            nyp2 = ny + 2
            nxp2 = nx + 2

            allocate(k2l2i(0:nx-1, 0:ny-1))
            allocate(k2l2(0:nx-1, 0:ny-1))

            allocate(filt(0:nz, 0:nx-1, 0:ny-1))

            allocate(a0(nx, ny))
            allocate(etdv(0:nz, nx, ny))
            allocate(htdv(0:nz, nx, ny))
            allocate(rkx(0:nx-1))
            allocate(hrkx(nx))
            allocate(rky(0:ny-1))
            allocate(hrky(ny))
            allocate(rkz(0:nz))
            allocate(rkzi(1:nz-1))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))
            allocate(ztrig(2 * nz))

            nxsub = nx / nsubs_tri

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(0) = zero
            do kx = 1, nwx-1
                rkx(kx)    = hrkx(2 * kx)
                rkx(nx-kx) = hrkx(2 * kx)
            enddo
            rkx(nwx) = hrkx(nx)
            rkxmax = hrkx(nx)

            !Define y wavenumbers:
            rky(0) = zero
            do ky = 1, nwy-1
                rky(ky)    = hrky(2 * ky)
                rky(ny-ky) = hrky(2 * ky)
            enddo
            rky(nwy) = hrky(ny)
            rkymax = hrky(ny)

            !Define z wavenumbers:
            rkz(0) = zero
            call init_deriv(nz, extent(3), rkz(1:nz))
            rkzi(1:nz-1) = one / rkz(1:nz-1)

            !Squared maximum total wavenumber:
            rksqmax = rkxmax ** 2 + rkymax ** 2

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
        end subroutine

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
            do ky = 0, ny-1
                do kx = 0, nx-1
                    call dst(1, nz, sfc(1:nz, kx, ky), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

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
            integer                         :: iz, kx, ky

            ! transform sf(1:nz-1, :, :) to semi-spectral space (sine transform) as the array sf:
            !$omp parallel workshare
            sftop = sf(nz, :, :)
            !$omp end parallel workshare

            !$omp parallel do collapse(2)
            do ky = 0, ny-1
                do kx = 0, nx-1
                    sf(nz, kx, ky) = zero
                    call dst(1, nz, sf(1:nz, kx, ky), ztrig, zfactors)
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

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs, ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: kx, dkx, kxc

            !Carry out differentiation by wavenumber multiplication:
            ds(:, 1, :) = zero
            do kx = 2, nx - nwx
                dkx = 2 * (kx - 1)
                kxc = nxp2 - kx
                ds(:, kx,  :) = -hrkx(dkx) * fs(:,kxc,:)
                ds(:, kxc, :) =  hrkx(dkx) * fs(:,kx ,:)
            enddo

            if (mod(nx, 2) .eq. 0) then
                kxc = nwx + 1
                ds(:, kxc, :) = zero
            endif
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dy
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffy(fs,ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            double precision              :: fac
            integer                       :: ky, kyc

            !Carry out differentiation by wavenumber multiplication:
            ds(:, :, 1) = zero

            do ky = 2, ny - nwy
                kyc = nyp2 - ky
                fac = hrky(2 * (ky - 1))
                ds(:, :, ky) = -fac * fs(:, :, kyc)
                ds(:, :, kyc) = fac * fs(:, : , ky)
            enddo

            if (mod(ny, 2) .eq. 0) then
                kyc = nwy + 1
                ds(:, :, kyc) = zero
            endif
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f using 2nd-order differencing.
        !Here fs = f, ds = df/dz.
        subroutine diffz(fs, ds)
            double precision, intent(in)  :: fs(0:nz, 0:ny-1, 0:nx-1)
            double precision, intent(out) :: ds(0:nz, 0:ny-1, 0:nx-1)
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

        end subroutine

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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        ! and returns the result as fs in spectral space (in x & y).
        ! Only FFTs over the x and y directions are performed.
        ! *** fp is destroyed upon exit ***
        subroutine fftxyp2s(fp, fs)
            double precision, intent(inout) :: fp(:, :, :)       !Physical
            double precision, intent(out)   :: fs(:, :, :)       !Spectral
            integer                         :: kx, iy, nzval, nxval, nyval

            nzval = size(fp, 1)
            nyval = size(fp, 2)
            nxval = size(fp, 3)

            ! Carry out a full x transform first:
            call forfft(nzval * nyval, nxval, fp, xtrig, xfactors)

            ! Transpose array:
            !$omp parallel do collapse(2) shared(fs, fp) private(kx, iy)
            do kx = 1, nxval
                do iy = 1, nyval
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full y transform on transposed array:
            call forfft(nzval * nxval, nyval, fs, ytrig, yfactors)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes an *inverse* 2D FFT (in x & y) of a 3D array fs in spectral
        ! space and returns the result as fp in physical space (in x & y).
        ! Only inverse FFTs over the x and y directions are performed.
        ! *** fs is destroyed upon exit ***
        subroutine fftxys2p(fs, fp)
            double precision, intent(inout):: fs(:, :, :)  !Spectral
            double precision, intent(out)  :: fp(:, :, :)  !Physical
            integer                        :: kx, iy, nzval, nxval, nyval

            nzval = size(fs, 1)
            nxval = size(fs, 2)
            nyval = size(fs, 3)

            ! Carry out a full inverse y transform first:
            call revfft(nzval * nxval, nyval, fs, ytrig, yfactors)

            ! Transpose array:
            !$omp parallel do collapse(2) shared(fs, fp) private(kx, iy)
            do kx = 1, nxval
                do iy = 1, nyval
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full inverse x transform:
            call revfft(nzval * nyval, nxval, fp, xtrig, xfactors)
        end subroutine

end module inversion_utils
