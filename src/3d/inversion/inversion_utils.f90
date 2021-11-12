module inversion_utils
    use constants
    use parameters, only : nx, ny, nz, dx, dxi, extent
    use stafft
    use sta2dfft
    implicit none

    private

    ! Ordering in physical space: z, y, x
    ! Ordering in spectral space: z, x, y

    ! Tridiagonal arrays for the horizontal vorticity components:
    double precision, allocatable :: etdh(:, :, :), htdh(:, :, :), ap(:, :), apb(:, :)

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    ! Tridiagonal arrays for the compact difference calculation of d/dz:
    double precision, allocatable :: etd0(:), htd0(:)
    double precision, allocatable :: etd1(:), htd1(:)

    ! Tridiagonal arrays for integrating in z:
    double precision, allocatable :: etda(:), htda(:)

    !Horizontal wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:)
    integer :: xfactors(5), yfactors(5)
    integer, parameter :: nsubs_tri = 8 !number of blocks for openmp
    integer :: nxsub

    !De-aliasing filter:
    double precision, allocatable :: filt(:, :)

    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq

    logical :: is_initialised = .false.

    public :: init_fft  &
            , diffx     &
            , diffy     &
            , diffz0    &
            , diffz1    &
            , lapinv0   &
            , lapinv1   &
            , vertint   &
            , fftxyp2s  &
            , fftxys2p  &
            , dz2       &
            , filt      &
            , hdzi      &
            , xfactors  &
            , yfactors  &
            , xtrig     &
            , ytrig

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision   :: a0(nx, ny), a0b(nx, ny), ksq(nx, ny)
            double precision   :: rkxmax, rkymax
            double precision   :: rksqmax, rkfsq
            integer            :: nwx, nwy
            integer            :: kx, ky, iz, isub, ib_sub, ie_sub

            if (is_initialised) then
                return
            endif

            is_initialised = .true.

            dz = dx(3)
            dzi = dxi(3)
            dz6  = f16 * dx(3)
            dz2  = f12 * dx(3)
            dz24 = f124 * dx(3)
            dzisq = dxi(3) ** 2
            hdzi = f12 * dxi(3)
            nwx = nx / 2
            nwy = ny / 2

            allocate(etdh(nz-1, nx, ny))
            allocate(htdh(nz-1, nx, ny))
            allocate(ap(nx, ny))
            allocate(apb(nx, ny))
            allocate(etdv(0:nz, nx, ny))
            allocate(htdv(0:nz, nx, ny))
            allocate(etd0(0:nz))
            allocate(htd0(0:nz))
            allocate(etd1(nz-1))
            allocate(htd1(nz-1))
            allocate(etda(nz))
            allocate(htda(nz))
            allocate(rkx(nx))
            allocate(hrkx(nx))
            allocate(rky(ny))
            allocate(hrky(ny))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))
            allocate(filt(nx, ny))

            nxsub = nx / nsubs_tri

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)

            !Define x wavenumbers:
            rkx(1) = zero
            do kx = 1, nwx-1
                rkx(kx+1)    = hrkx(2 * kx)
                rkx(nx+1-kx) = hrkx(2 * kx)
            enddo
            rkx(nwx+1) = hrkx(nx)
            rkxmax = hrkx(nx)

            !Define y wavenumbers:
            rky(1) = zero
            do ky = 1, nwy-1
                rky(ky+1)    = hrky(2 * ky)
                rky(ny+1-ky) = hrky(2 * ky)
            enddo
            rky(nwy+1) = hrky(ny)
            rkymax = hrky(ny)

            !Squared maximum total wavenumber:
            rksqmax = rkxmax ** 2 + rkymax ** 2

            !Squared wavenumber array (used in tridiagonal solve):
            do ky = 1, ny
                do kx = 1, nx
                    ksq(kx, ky) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            !--------------------------------------------------------------------
            ! Define de-aliasing filter:
            rkfsq = two * rksqmax / nine
            ! rkfsq: the square of the filter wavenumber (generic 2/3 rule)
            do ky = 1, ny
                do kx = 1, nx
                    if (ksq(kx, ky) .gt. rkfsq) then
                        filt(kx, ky) = zero
                    else
                        filt(kx, ky) = one
                    endif
                enddo
            enddo

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0 = -two * dzisq - f56 * ksq
            a0b = -dzisq - f13 * ksq
            ap = dzisq - f112 * ksq
            apb = dzisq - f16 * ksq

            !-----------------------------------------------------------------------
            ! Tridiagonal arrays for the horizontal vorticity components:
            htdh(1, :, :) = one / a0
            etdh(1, :, :) = -ap * htdh(1, :, :)
            !$omp parallel shared(a0, ap, etdh, htdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-2
                    htdh(iz, ib_sub:ie_sub, :) = one / (a0(ib_sub:ie_sub, :) &
                                               + ap(ib_sub:ie_sub, :) * etdh(iz-1, ib_sub:ie_sub, :))
                    etdh(iz, ib_sub:ie_sub, :) = -ap(ib_sub:ie_sub, :) * htdh(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
            htdh(nz-1, :, :) = one / (a0 + ap * etdh(nz-2, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdh(:, 1, 1) = zero
            etdh(:, 1, 1) = zero

            ! Tridiagonal arrays for the vertical vorticity component:
            htdv(0, :, :) = one / a0b
            etdv(0, :, :) = -apb * htdv(0, :, :)
            !$omp parallel shared(a0, ap, etdv, htdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    htdv(iz, ib_sub:ie_sub, :) = one / (a0(ib_sub:ie_sub, :) &
                                               + ap(ib_sub:ie_sub, :) * etdv(iz-1, ib_sub:ie_sub, :))
                    etdv(iz, ib_sub:ie_sub, :) = -ap(ib_sub:ie_sub, :) * htdv(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            etdv(nz-1, 1, 1) = zero

            htdv(nz, :, :) = one / (a0b + apb * etdv(nz-1, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdv(:, 1, 1) = zero
            etdv(:, 1, 1) = zero

            ! Tridiagonal arrays for the compact difference calculation of d/dz
            ! for fields f for which f  =  0 at the boundaries:
            htd0(0) = one / f23
            etd0(0) = -f13 * htd0(0)
            do iz = 1, nz-1
                htd0(iz) = one / (f23 + f16 * etd0(iz-1))
                etd0(iz) = -f16 * htd0(iz)
            enddo
            htd0(nz) = one / (f23 + f13 * etd0(nz-1))

            ! Tridiagonal arrays for the compact difference calculation of d /dz
            ! for fields f for which df/dz = 0 at the boundaries:
            htd1(1) = one / f23
            etd1(1) = -f16 * htd1(1)
            do iz = 2, nz-2
                htd1(iz) = one / (f23 + f16 * etd1(iz-1))
                etd1(iz) = -f16 * htd1(iz)
            enddo
            htd1(nz-1) = one / (f23 + f16 * etd1(nz-2))

            ! Tridiagonal arrays used for integrating in z (see vertint):
            htda(1) = one / (one + f16)
            etda(1) = -f16 * htda(1)
            do iz = 2, nz-1
                htda(iz) = one / (one + f16 * etda(iz-1))
                etda(iz) = -f16 * htda(iz)
            enddo
            htda(nz) = one / (one + f16 + f16 * etda(nz-2))
        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs,ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: nwx, nxp2, kx, dkx, kxc

            nwx = nx / 2
            nxp2 = nx + 2

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
            integer                       :: nwy, nyp2, ky, kyc

            !Could be pre-defined constants:
            nwy = ny / 2
            nyp2 = ny + 2

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

        !Calculates df/dz for a field f which has f = 0 at the boundaries
        !using 4th-order compact differencing.  Here fs = f, ds = df/dz,
        !lapfsbot = Laplace(f) at z_min and lapfstop = Laplace(f) at z_max,
        !both given.  *** All quantities must be in semi-spectral space ***
        subroutine diffz0(fs, ds, lapfsbot, lapfstop)
            double precision, intent(in)  :: fs(0:nz, nx, ny), &
                                             lapfsbot(nx, ny), &
                                             lapfstop(nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: iz, isub, ib_sub, ie_sub

            ds(0, :, :) = fs(1, :, :) * dzi - dz6 * lapfsbot
            ds(1, :, :) = fs(2, :, :) * hdzi

            !$omp parallel shared(ds, fs, hdzi, nz) private(iz) default(none)
            !$omp do
            do iz = 2, nz-2
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz-1, :, :) =  - fs(nz-2, :, :) * hdzi
            ds(nz, :, :) = f16 * dz * lapfstop - fs(nz-1, :, :) * dzi

            ds(0, :, :) = ds(0, :, :) * htd0(0)

            !$omp parallel shared(ds, htd0, nz, nxsub) private(isub, ib_sub, ie_sub,iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub+1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    ds(iz, ib_sub:ie_sub, :) = (ds(iz, ib_sub:ie_sub, :) &
                                             - f16 * ds(iz-1, ib_sub:ie_sub, :)) * htd0(iz)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz, :, :) = (ds(nz, :, :) - f13 * ds(nz-1, :, :)) * htd0(nz)

            !$omp parallel shared(ds, etd0, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-1, 0, -1
                    ds(iz, ib_sub:ie_sub, :) = etd0(iz) * ds(iz+1, ib_sub:ie_sub, :) &
                                             + ds(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f which has df/dz = 0 at the boundaries
        !using 4th-order compact differencing.  Here fs = f and ds = df/dz.
        subroutine diffz1(fs, ds)
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: iz, isub, ib_sub, ie_sub

            !$omp parallel shared(ds, fs, nz, hdzi) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

            ds(0, :, :) = zero
            ds(1, :, :) = ds(1, :, :) * htd1(1)

            !$omp parallel shared(ds, htd1, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    ds(iz, ib_sub:ie_sub, :) = (ds(iz, ib_sub:ie_sub, :) &
                                             - f16 * ds(iz-1, ib_sub:ie_sub, :)) * htd1(iz)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz, :, :) = zero

            !$omp parallel shared(ds, etd1, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-2, 1, -1
                    ds(iz, ib_sub:ie_sub, :) = etd1(iz) * ds(iz+1, ib_sub:ie_sub, :) &
                                             + ds(iz, ib_sub:ie_sub, :)
                enddo
            end do
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here fs = 0 on the z boundaries.
        !Uses 4th-order compact differencing
        !*** Overwrites fs ***
        subroutine lapinv0(fs)
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            double precision                :: rs(nz-1, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            !$omp parallel shared(rs, fs, nz) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                rs(iz, :, :) = f112 * (fs(iz-1, :, :) + fs(iz+1, :, :)) + f56 * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            fs(0, :, :) = zero
            fs(1, :, :) = rs(1, :, :) * htdh(1, :, :)

            !$omp parallel shared(rs, fs, ap, htdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    fs(iz, ib_sub:ie_sub, :) = (rs(iz, ib_sub:ie_sub, :) &
                                             - ap(ib_sub:ie_sub, :) * fs(iz-1, ib_sub:ie_sub, :)) &
                                             * htdh(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = zero

            !$omp parallel shared(fs, etdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz)  default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-2, 1, -1
                    fs(iz, ib_sub:ie_sub, :) = etdh(iz, ib_sub:ie_sub, :) &
                                             * fs(iz+1, ib_sub:ie_sub, :) + fs(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here dfs/dz = 0 on the z boundaries.
        !Uses 4th-order compact differencing
        !*** Overwrites fs ***
        subroutine lapinv1(fs)
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            double precision                :: rs(0:nz, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            rs(0, :, :) = f13 * fs(0, :, :) + f16 * fs(1, :, :)
            !$omp parallel shared(rs, fs) private(iz)
            !$omp do
            do iz = 1, nz-1
                rs(iz, :, :) = f112 * (fs(iz-1, :, :) + fs(iz+1, :, :)) + f56 * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel
            rs(nz, :, :) = f13 * fs(nz, :, :) + f16 * fs(nz-1, :, :)

            fs(0, :, :) = rs(0, :, :) * htdv(0, :, :)

            !$omp parallel shared(rs, fs, ap, htdv, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    fs(iz, ib_sub:ie_sub, :) = (rs(iz, ib_sub:ie_sub, :) &
                                             - ap(ib_sub:ie_sub, :) * fs(iz-1, ib_sub:ie_sub, :)) &
                                             * htdv(iz, ib_sub:ie_sub, :)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = (rs(nz, :, :) - apb * fs(nz-1, :, :)) * htdv(nz, :, :)

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

        !Finds f by integrating df/dz = d, ensuring f = 0 at the boundaries
        !using 4th-order compact differencing.  Here ds = df/dz and fs = f.
        subroutine vertint(ds, fs)
            double precision, intent(in)  :: ds(0:nz)
            double precision, intent(out) :: fs(0:nz)
            double precision              :: es(nz), esum
            integer                       :: iz

            !-------------------------------------------
            !First interpolate ds to a half grid as es:
            do iz = 1, nz
                es(iz) = f23 * (ds(iz-1) + ds(iz))
            enddo

            es(1) = es(1) * htda(1)
            do iz = 2, nz
                es(iz) = (es(iz) - f16 * es(iz-1)) * htda(iz)
            enddo

            do iz = nz-1, 1, -1
                es(iz) = etda(iz) * es(iz+1) + es(iz)
            enddo

            !-------------------------------------------
            !Next adjust es to ensure f(nz) = 0:
            esum = (f1112 * (es(1) + es(nz)) + sum(es(2:nz-1))) / (dble(nz) - f16)
            es = es - esum

            !Integrate:
            fs(0) = zero
            fs(1) = dz24 * (23.0d0 * es(1) + es(2))
            do iz = 2, nz-1
                fs(iz) = fs(iz-1) + dz24 * (es(iz-1) + 22.d0 * es(iz) + es(iz+1))
            enddo
            fs(nz) = zero
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        ! and returns the result as fs in spectral space (in x & y).
        ! Only FFTs over the x and y directions are performed.
        ! *** fp is destroyed upon exit ***
        subroutine fftxyp2s(fp, fs)
            double precision, intent(inout) :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out)   :: fs(0:nz, nx, ny)  !Spectral
            integer                         :: kx, iy

            ! Carry out a full x transform first:
            call forfft((nz+1) * ny, nx, fp, xtrig, xfactors)

            ! Transpose array:
            !$omp parallel do
            do kx = 1, nx
                do iy = 1, ny
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full y transform on transposed array:
            call forfft((nz+1) * nx, ny, fs, ytrig, yfactors)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes an *inverse* 2D FFT (in x & y) of a 3D array fs in spectral
        ! space and returns the result as fp in physical space (in x & y).
        ! Only inverse FFTs over the x and y directions are performed.
        ! *** fs is destroyed upon exit ***
        subroutine fftxys2p(fs,fp)
            double precision, intent(inout):: fs(0:nz, nx, ny)  !Spectral
            double precision, intent(out)  :: fp(0:nz, ny, nx)  !Physical
            integer                        :: kx, iy

            ! Carry out a full inverse y transform first:
            call revfft((nz+1) * nx, ny, fs, ytrig, yfactors)

            ! Transpose array:
            !$omp parallel do
            do kx = 1, nx
                do iy = 1, ny
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo
            !$omp end parallel do

            ! Carry out a full inverse x transform:
            call revfft((nz+1) * ny, nx, fp, xtrig, xfactors)
        end subroutine

end module inversion_utils
