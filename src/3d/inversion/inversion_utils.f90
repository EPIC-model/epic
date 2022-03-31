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
    double precision, allocatable :: etdh(:, :, :), htdh(:, :, :)

    ! Tridiagonal arrays for the vertical vorticity component:
    double precision, allocatable :: etdv(:, :, :), htdv(:, :, :)

    !Horizontal wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:)
    integer :: xfactors(5), yfactors(5)
    integer, parameter :: nsubs_tri = 8 !number of blocks for openmp
    integer :: nxsub




    double precision :: dz, dzi, dz2, dz6, dz24, hdzi, dzisq, ap, a0b
    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_initialised = .false.

    public :: init_fft  &
            , diffx     &
            , diffy     &
            , diffz     &
            , lapinv0   &
            , lapinv1   &
            , vertint   &
            , fftxyp2s  &
            , fftxys2p  &
            , dz2       &
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
            double precision, allocatable  :: a0(:, :), ksq(:, :)
            double precision               :: rkxmax, rkymax
            double precision               :: rksqmax
            integer                        :: kx, ky, iz, isub, ib_sub, ie_sub

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
            nyp2 = ny + 2
            nxp2 = nx + 2

            allocate(a0(nx, ny))
            allocate(ksq(nx, ny))

            allocate(etdh(nz-1, nx, ny))
            allocate(htdh(nz-1, nx, ny))
            allocate(etdv(0:nz, nx, ny))
            allocate(htdv(0:nz, nx, ny))
            allocate(rkx(nx))
            allocate(hrkx(nx))
            allocate(rky(ny))
            allocate(hrky(ny))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))

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

            !-----------------------------------------------------------------------
            ! Fixed coefficients used in the tridiagonal problems:
            a0 = -two * dzisq - ksq
            a0b = -dzisq
            ap = dzisq

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
                                               + ap * etdh(iz-1, ib_sub:ie_sub, :))
                    etdh(iz, ib_sub:ie_sub, :) = -ap * htdh(iz, ib_sub:ie_sub, :)
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
            etdv(0, :, :) = -ap * htdv(0, :, :)
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

            htdv(nz, :, :) = one / (a0b + ap * etdv(nz-1, :, :))
            ! Remove horizontally-averaged part (done separately):
            htdv(:, 1, 1) = zero
            etdv(:, 1, 1) = zero

            deallocate(a0)
            deallocate(ksq)
        end subroutine


        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs,ds)
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
            double precision, intent(in)  :: fs(0:nz, nx, ny)
            double precision, intent(out) :: ds(0:nz, nx, ny)
            integer                       :: iz

            ! forward and backward differencing for boundary cells
            ! iz = 0:  (fs(1) - fs(0)) / dz
            ! iz = nz: (fs(nz) - fs(nz-1)) / dz
            ds(0,  :, :) =   dzi * (fs(1,    :, :) - fs(0,    :, :))
            ds(nz, :, :) = - dzi * (fs(nz,   :, :) - fs(nz-1, :, :))

            ! central differencing for interior cells
            !$omp parallel shared(ds, fs, hdzi, nz) private(iz) default(none)
            !$omp do
            do iz = 1, nz-1
                ds(iz, :, :) = (fs(iz+1, :, :) - fs(iz-1, :, :)) * hdzi
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Inverts Laplace's operator on fs in semi-spectral space.
        !Here fs = 0 on the z boundaries.
        !Uses 2nd-order differencing
        !*** Overwrites fs ***
        subroutine lapinv0(fs)
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            fs(0, :, :) = zero
            fs(1, :, :) = fs(1, :, :) * htdh(1, :, :)

            !$omp parallel shared(fs, ap, htdh, nz, nxsub) private(isub, ib_sub, ie_sub, iz) default(none)
            !$omp do
            do isub = 0, nsubs_tri-1
                ib_sub = isub * nxsub + 1
                ie_sub = (isub + 1) * nxsub
                do iz = 2, nz-1
                    fs(iz, ib_sub:ie_sub, :) = (fs(iz, ib_sub:ie_sub, :) &
                                             - ap * fs(iz-1, ib_sub:ie_sub, :)) &
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
        !Uses 2nd-order differencing
        !*** Overwrites fs ***
        subroutine lapinv1(fs)
            double precision, intent(inout) :: fs(0:nz, nx, ny)
            double precision                :: rs(0:nz, nx, ny)
            integer                         :: iz, isub, ib_sub, ie_sub

            fs(0, :, :) = fs(0, :, :) * htdv(0, :, :)
            rs = fs

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

            fs(nz, :, :) = (rs(nz, :, :) - ap * fs(nz-1, :, :)) * htdv(nz, :, :)

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
        !using finite differencing.  Here ds = df/dz and fs = f.
        subroutine vertint(ds, fs)
            double precision, intent(in)  :: ds(0:nz)
            double precision, intent(out) :: fs(0:nz)
            double precision              :: c
            integer                       :: iz

            ! set lower boundary value
            fs(0)  = zero

            !$omp parallel private(iz)
            !$omp do
            do iz = 1, nz
                fs(iz) = fs(iz-1) + dz2 * (ds(i) - ds(i-1))
            enddo
            !$omp end do
            !$omp end parallel

            ! shift to adjust f(nz) to be zero
            c = f(nz) / dble(nz)

            !$omp parallel private(iz)
            !$omp do
            do iz = 1, nz
                fs(iz) = fs(iz) - c * dble(iz)
            enddo
            !$omp end do
            !$omp end parallel

            ! set upper boundary value
            fs(nz)  = zero

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
            !$omp parallel do
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
            !$omp parallel do
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
