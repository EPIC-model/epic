module fft
    use constants
    use parameters
    use stafft
    use sta2dfft
    implicit none

    ! Tridiagonal arrays for the horizontal vorticity components:
    double precision :: etdh(nz-1, ny, nx), htdh(nz-1, ny, nx), ap(ny, nx), apb(ny, nx)

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision :: etdv(0:nz, ny, nx), htdv(0:nz, ny, nx)

    ! Tridiagonal arrays for the compact difference calculation of d/dz:
    double precision :: etd0(0:nz), htd0(0:nz)
    double precision :: etd1(nz-1), htd1(nz-1)

    ! Tridiagonal arrays for integrating in z:
    double precision :: etda(nz), htda(nz)

    !Horizontal wavenumbers:
    double precision :: rkx(nx), hrkx(nx), rky(ny), hrky(ny)

    !Quantities needed in FFTs:
    double precision :: xtrig(2 * nx), ytrig(2 * ny)
    integer :: xfactors(5), yfactors(5)
    integer, parameter :: nsubs_tri = 8 !number of blocks for openmp
    integer, parameter :: nxsub = nx / nsubs_tri !number of x cells per block

    !De-aliasing filter:
    double precision :: filt(nx, ny)

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Initialises this module (FFTs, x & y wavenumbers, tri-diagonal
        !coefficients, etc).
        subroutine init_fft
            double precision   :: a0(ny, nx), a0b(ny, nx), ksq(ny, nx)
            double precision   :: rkxmax, rkymax
            double precision   :: rksqmax, rkfsq
            integer, parameter :: nwx = nx / 2,nwy = ny / 2
            integer            :: kx, ky, iz, isub, ib_sub, ie_sub

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, ellx, elly, xfactors, yfactors, xtrig, ytrig, hrkx, hrky)

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
            do kx = 1, nx
                do ky = 1, ny
                    ksq(ky, kx) = rkx(kx) ** 2 + rky(ky) ** 2
                enddo
            enddo

            !--------------------------------------------------------------------
            ! Define de-aliasing filter:
            rkfsq = two * rksqmax / nine
            ! rkfsq: the square of the filter wavenumber (generic 2/3 rule)
            do kx = 1, nx
                do ky = 1, ny
                    if (ksq(ky, kx) .gt. rkfsq) then
                        filt(ky, kx) = zero
                    else
                        filt(ky, kx) = one
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
            !$omp parallel shared(a0, ap, etdh, htdh) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-2
                    htdh(iz, :, ib_sub:ie_sub) = one / (a0(:, ib_sub:ie_sub) &
                                               + ap(:, ib_sub:ie_sub) * etdh(iz-1, :, ib_sub:ie_sub))
                    etdh(iz, :, ib_sub:ie_sub) = -ap(:, ib_sub:ie_sub) * htdh(iz, :, ib_sub:ie_sub)
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
            !$omp parallel shared(a0, ap, etdv, htdv) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    htdv(iz, :, ib_sub:ie_sub) = one / (a0(:, ib_sub:ie_sub) &
                                               + ap(:, ib_sub:ie_sub) * etdv(iz-1, :, ib_sub:ie_sub))
                    etdv(iz, :, ib_sub:ie_sub) = -ap(:, ib_sub:ie_sub) * htdv(iz, :, ib_sub:ie_sub)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

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

        !Given fs in spectral space (at least in x & y), this returns dfs/dx
        !(partial derivative).  The result is returned in ds, again
        !spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs, ds)
            double precision, intent(in)  :: fs(0:nz, ny, nx)
            double precision, intent(out) :: ds(0:nz, ny, nx)
            integer                       :: iz

            !$omp parallel shared(ds,fs) private(iz) firstprivate(hrkx) default(none)
            !$omp do
            do iz = 0, nz
                call xderiv(nx, ny, hrkx, fs(iz, 1, 1), ds(iz, 1, 1))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Given fs in spectral space (at least in x & y), this returns dfs/dy
        !(partial derivative).  The result is returned in ds, again
        !spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffy(fs, ds)
            double precision, intent(in)  :: fs(0:nz, ny, nx)
            double precision, intent(out) :: ds(0:nz, ny, nx)
            integer                       :: iz

            !$omp parallel shared(ds,fs) private(iz) firstprivate(hrky) default(none)
            !$omp do
            do iz = 0, nz
                call yderiv(nx, ny, hrky, fs(iz, 1, 1), ds(iz, 1, 1))
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f which has f = 0 at the boundaries
        !using 4th-order compact differencing.  Here fs = f, ds = df/dz,
        !lapfsbot = Laplace(f) at z_min and lapfstop = Laplace(f) at z_max,
        !both given.  *** All quantities must be in semi-spectral space ***
        subroutine diffz0(fs, ds, lapfsbot, lapfstop)
            double precision, intent(in)  :: fs(0:nz, ny, nx), &
                                             lapfsbot(ny, nx), &
                                             lapfstop(ny, nx)
            double precision, intent(out) :: ds(0:nz, ny, nx)
            integer                       :: iz, isub, ib_sub, ie_sub

            ds(0, :, :) = fs(1, :, :) * dzi - dz6 * lapfsbot
            ds(1, :, :) = fs(2, :, :) * hdzi

            !$omp parallel shared(ds,fs) private(iz) default(none)
            !$omp do
            do iz = 2, nz-2
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz-1, :, :) =  -fs(nz-2:, :) * hdzi
            ds(nz, :, :) = f16 * dz * lapfstop - fs(nz-1, :, :) * dzi

            ds(0, :, :) = ds(0, :, :) * htd0(0)

            !$omp parallel shared(ds,htd0) private(isub,ib_sub,ie_sub,iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub+1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    ds(iz, :, ib_sub:ie_sub) = (ds(iz, :, ib_sub:ie_sub) &
                                             - f16 * ds(iz-1, :, ib_sub:ie_sub)) * htd0(iz)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz, :, :) = (ds(nz, :, :) - f13 * ds(nz-1, :, :)) * htd0(nz)

            !$omp parallel shared(ds,etd0) private(isub,ib_sub,ie_sub,iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-1, 0, -1
                    ds(iz, :, ib_sub:ie_sub) = etd0(iz) * ds(iz+1, :, ib_sub:ie_sub) &
                                             + ds(iz, :, ib_sub:ie_sub)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Calculates df/dz for a field f which has df/dz = 0 at the boundaries
        !using 4th-order compact differencing.  Here fs = f and ds = df/dz.
        subroutine diffz1(fs, ds)
            double precision, intent(in)  :: fs(0:nz, ny, nx)
            double precision, intent(out) :: ds(0:nz, ny, nx)
            integer                       :: iz, isub, ib_sub, ie_sub

            !$omp parallel shared(ds, fs) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

            ds(0, :, :) = zero
            ds(1, :, :) = ds(1, :, :) * htd1(1)

            !$omp parallel shared(ds, htd1) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    ds(iz, :, ib_sub:ie_sub) = (ds(iz, :, ib_sub:ie_sub)
                                             - f16 * ds(iz-1, :, ib_sub:ie_sub)) * htd1(iz)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            ds(nz, :, :) = zero

            !$omp parallel shared(ds, etd1) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-2, 1, -1
                    ds(iz, :, ib_sub:ie_sub) = etd1(iz) * ds(iz+1, :, ib_sub:ie_sub) &
                                             + ds(iz, :, ib_sub:ie_sub)
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
            double precision, intent(inout) :: fs(0:nz, ny, nx)
            double precision                :: rs(nz-1, ny, nx)
            integer                         :: iz, isub, ib_sub, ie_sub

            !$omp parallel shared(rs, fs) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                rs(iz, :, :) = f112 * (fs(iz-1, :, :) + fs(iz+1, :, :)) + f56 * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            fs(0, :, :) = zero
            fs(1, :, :) = rs(1, :, :) * htdh(1, :, :)

            !$omp parallel shared(rs, fs, ap, htdh) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    fs(iz, :, ib_sub:ie_sub) = (rs(iz, :, ib_sub:ie_sub) &
                                             - ap(:, ib_sub:ie_sub) * fs(iz-1, :, ib_sub:ie_sub)) &
                                             * htdh(iz, :, ib_sub:ie_sub)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = zero

            !$omp parallel shared(fs, etdh) private(isub, ib_sub, ie_sub, iz)  default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-2, 1, -1
                    fs(iz, :, ib_sub:ie_sub) = etdh(iz, :, ib_sub:ie_sub) &
                                             * fs(iz+1, :, ib_sub:ie_sub) + fs(iz, :, ib_sub:ie_sub)
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
            double precision, intent(inout) :: fs(0:nz, ny, nx)
            double precision                :: rs(0:nz, ny, nx)
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

            !$omp parallel shared(rs, fs, ap, htdv) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 1, nz-1
                    fs(iz, :, ib_sub:ie_sub) = (rs(iz, :, ib_sub:ie_sub) &
                                             - ap(:, ib_sub:ie_sub) * fs(iz-1, :, ib_sub:ie_sub)) &
                                             * htdv(iz, :, ib_sub:ie_sub)
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            fs(nz, :, :) = (rs(nz, :, :) - apb * fs(nz-1, :, :)) * htdv(nz, :, :)

            !$omp parallel shared(fs, etdv) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = nz-1, 0, -1
                    fs(iz, :, ib_sub:ie_sub) = etdv(iz, :, ib_sub:ie_sub) * fs(iz+1, :, ib_sub:ie_sub) &
                                             + fs(iz, :, ib_sub:ie_sub)
                enddo
            enddo
            !$omp end do
            !$omp end parallel
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
                fs(iz) = fs(iz-1) + dz24 * (es(iz-1) + 22.0d0 * es(iz) + es(iz+1))
            enddo
            fs(nz) = zero
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        !and returns the result as fs in spectral space (in x & y).
        !Only FFTs over the x and y directions are performed.
        subroutine fftxyp2s(fp, fs)
            double precision, intent(in)  :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out) :: fs(0:nz, nx, ny)  !Spectral
            double precision              :: fpz(ny, nx)  !Physical
            double precision              :: fsz(nx, ny)  !Spectral
            integer                       :: iz

            !$omp parallel shared(fp, fs) private(iz, fpz, fsz) &
            !$omp& firstprivate(xfactors, yfactors, xtrig, ytrig) default(none)
            !$omp do
            do iz = 0, nz
                fpz = fp(iz, :, :)
                call ptospc(nx, ny, fpz, fsz, xfactors, yfactors, xtrig, ytrig)
                fs(iz, :, :) = fsz
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes an *inverse* 2D FFT (in x & y) of a 3D array fs in spectral
        !space and returns the result as fp in physical space (in x & y).
        !Only inverse FFTs over the x and y directions are performed.
        subroutine fftxys2p(fs, fp)
            double precision, intent(out) :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(in)  :: fs(0:nz, nx, ny)  !Spectral
            double precision              :: fpz(ny, nx)  !Physical
            double precision              :: fsz(nx, ny)  !Spectral
            integer                       :: iz

            !$omp parallel shared(fp, fs) private(iz, fpz, fsz) &
            !$omp& firstprivate(xfactors, yfactors, xtrig, ytrig) default(none)
            !$omp do
            do iz = 0, nz
                fsz = fs(iz, :, :)
                call spctop(nx, ny, fsz, fpz, xfactors, yfactors, xtrig, ytrig)
                fp(iz, :, :) = fpz
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Given the vorticity vector (ox, oy, oz) in physical space, this
        !returns the associated velocity field (uu, vv, ww) and a spectral
        !copy (us, vs, ws),  which is additionally filtered.  Note: the
        !vorticity is modified to be solenoidal and spectrally filtered.
        subroutine vor2vel(ox, oy, oz,  uu, vv, ww,  us, vs, ws)
            double precision, intent(in)  :: ox(0:nz, ny, nx), oy(0:nz, ny, nx), oz(0:nz, ny, nx)
            double precision, intent(out) :: uu(0:nz, ny, nx), vv(0:nz, ny, nx), ww(0:nz, ny, nx)
            double precision, intent(out) :: us(0:nz, ny, nx), vs(0:nz, ny, nx), ws(0:nz, ny, nx)
            double precision              :: as(0:nz, ny, nx), bs(0:nz, ny, nx), cs(0:nz, ny, nx)
            double precision              :: ds(0:nz, ny, nx), es(0:nz, ny, nx), fs(0:nz, ny, nx)
            double precision              :: astop(ny, nx), bstop(ny, nx)
            double precision              :: asbot(ny, nx), bsbot(ny, nx)
            double precision              :: ubar(0:nz), vbar(0:nz)
            integer                       :: iz

            !------------------------------------------------------------------
            !Convert vorticity to spectral space as (as, bs, cs):
            call fftxyp2s(ox, as)
            call fftxyp2s(oy, bs)
            call fftxyp2s(oz, cs)

            !Add -grad(lambda) where Laplace(lambda) = div(ox, oy, oz) to
            !enforce the solenoidal condition on the vorticity field:
            call diffx(as, ds)
            call diffy(bs, es)

            !For the vertical parcel vorticity, use 2nd-order accurate
            !differencing with extrapolation at the boundaries:
            !$omp parallel shared(fs, cs) private(iz)
            !$omp do
            do iz = 1, nz-1
                fs(iz, :, :) = hdzi * (cs(iz+1, :, :) - cs(iz-1, :, :))
            enddo
            !$omp end do
            !$omp end parallel

            fs(0, :, :) = hdzi * (four * cs(1, :, :) - cs(2, :, :) - three * cs(0, :, :))
            fs(nz, :, :) = hdzi * (cs(nz-2, :, :) + three * cs(nz, :, :) - four * cs(nz-1, :, :))

            !Form div(ox, oy, oz):
            !$omp parallel
            !$omp workshare
            fs = ds + es + fs
            !$omp end workshare
            !$omp end parallel

            !Remove horizontally-averaged part (plays no role):
            fs(:, 1, 1) = zero

            !Invert Lap(lambda) = div(ox, oy, oz) assuming dlambda/dz = 0 at the
            !boundaries (store solution lambda in fs):
            call lapinv1(fs)

            !Filter lambda:
            !$omp parallel shared(fs, filt) private(iz)  default(none)
            !$omp do
            do iz = 0, nz
                fs(iz, :, :) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Subtract grad(lambda) to enforce div(ox, oy, oz) = 0:
            call diffx(fs, ds)
            !$omp parallel
            !$omp workshare
            as = as - ds
            !$omp end workshare
            !$omp end parallel

            call diffy(fs, ds)
            !$omp parallel
            !$omp workshare
            bs = bs - ds
            !$omp end workshare
            !$omp end parallel

            call diffz1(fs, ds)
            !$omp parallel
            !$omp workshare
            cs = cs - ds
            !$omp end workshare
            !$omp end parallel
            !Ensure horizontal average of vertical vorticity is zero:
            cs(:, 1, 1) = zero

            !Compute spectrally filtered vorticity in physical space:
            !$omp parallel shared(ds, es, fs, as, bs, cs, filt) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                ds(iz, :, :) = filt * as(iz, :, :)
                es(iz, :, :) = filt * bs(iz, :, :)
                fs(iz, :, :) = filt * cs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Save boundary values of x and y vorticity for z derivatives of A & B:
            asbot = as(0, :, :)
            bsbot = bs(0, :, :)
            astop = as(nz, :, :)
            bstop = bs(nz, :, :)

            !Return corrected vorticity to physical space in (ox, oy, oz):
            call fftxys2p(ds, ox)
            call fftxys2p(es, oy)
            call fftxys2p(fs, oz)

            !Define horizontally-averaged flow by integrating horizontal vorticity:
            ubar(0) = zero
            vbar(0) = zero
            do iz = 0, nz-1
                ubar(iz+1) = ubar(iz) + dz2 * (es(iz, 1, 1) + es(iz+1, 1, 1))
                vbar(iz+1) = vbar(iz) - dz2 * (ds(iz, 1, 1) + ds(iz+1, 1, 1))
            enddo

            !-----------------------------------------------------------------
            !Invert vorticity to find vector potential (A, B, C) -> (as, bs, cs):
            call lapinv0(as)
            call lapinv0(bs)
            call lapinv1(cs)

            !------------------------------------------------------------
            !Compute x velocity component, u = B_z - C_y:
            call diffy(cs, ds)
            call diffz0(bs, es, bsbot, bstop)
            !bsbot & bstop contain spectral y vorticity component at z_min and z_max
            !$omp parallel
            !$omp workshare
            fs = es - ds
            !$omp end workshare
            !$omp end parallel
            !Add horizontally-averaged flow:
            fs(:, 1, 1) = ubar
            !$omp parallel shared(us, fs, filt) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                us(iz, :, :) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get u in physical space:
            call fftxys2p(fs, uu)

            !------------------------------------------------------------
            !Compute y velocity component, v = C_x - A_z:
            call diffx(cs, ds)
            call diffz0(as, es, asbot, astop)
            !asbot & astop contain spectral x vorticity component at z_min and z_max
            !$omp parallel
            !$omp workshare
            fs = ds - es
            !$omp end workshare
            !$omp end parallel
            !Add horizontally-averaged flow:
            fs(:, 1, 1) = vbar
            !$omp parallel shared(vs, fs, filt) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                vs(iz, :, :) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get v in physical space:
            call fftxys2p(fs, vv)

            !------------------------------------------------------------
            !Compute z velocity component, w = A_y - B_x:
            call diffx(bs, ds)
            call diffy(as, es)
            !$omp parallel
            !$omp workshare
            fs = es - ds
            !$omp end workshare
            !$omp end parallel

            !$omp parallel shared(ws, fs, filt) private(iz) default(none)
            !$omp do
            do iz = 0, nz
                ws(iz, :, :) = filt * fs(iz, :, :)
            enddo
            !$omp end do
            !$omp end parallel

            !Get w in physical space:
            call fftxys2p(fs, ww)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a divergent flow field (ud, vd, wd) = grad(phi) where
        !Lap(phi) = div (given).
        subroutine diverge(div,  ud, vd, wd)
            double precision, intent(in)  :: div(0:nz, ny, nx)
            double precision, intent(out) :: ud(0:nz, ny, nx), vd(0:nz, ny, nx), wd(0:nz, ny, nx)
            double precision              :: ds(0:nz, ny, nx)
            double precision              :: us(0:nz, ny, nx), vs(0:nz, ny, nx), ws(0:nz, ny, nx)
            double precision              :: ads(0:nz), wbar(0:nz)
            integer                       :: iz

            !------------------------------------------------------------------
            !Convert div to spectral space (in x & y) as ds:
            call fftxyp2s(div, ds)

            ! Compute the x & y-independent part of ds by integration:
            do iz = 0, nz
                ads(iz) = ds(iz, 1, 1)
            enddo
            call vertint(ads, wbar)

            ! Invert Laplace's operator semi-spectrally with compact differences:
            call lapinv1(ds)

            ! Compute x derivative spectrally:
            call diffx(ds, us)

            ! Reverse FFT to define x velocity component ud:
            call fftxys2p(us, ud)

            ! Compute y derivative spectrally:
            call diffy(ds, vs)

            ! Reverse FFT to define y velocity component vd:
            call fftxys2p(vs, vd)

            ! Compute z derivative by compact differences:
            call diffz1(ds, ws)

            ! Add on the x and y-independent part of wd:
            do iz = 0, nz
                ws(iz, 1, 1) = ws(iz, 1, 1) + wbar(iz)
            enddo

            ! Reverse FFT to define z velocity component wd:
            call fftxys2p(ws, wd)
        end subroutine
end module fft
