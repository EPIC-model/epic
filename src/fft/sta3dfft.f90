module sta3dfft
    use fft_pencil, only : perform_fftxyp2s, perform_fftxys2p
    use mpi_layout
    use stafft, only : dct, dst
!     use sta2dfft
!     use deriv1d, only : init_deriv
    implicit none

    ! Wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:), rkzi(:)

    !Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer :: xfactors(5), yfactors(5), zfactors(5)

    integer :: nwx, nwy, nxp2, nyp2

    logical :: is_fft_initialised = .false.

    private :: is_fft_initialised

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_fft(nx, ny, nz, lx, ly, lz)
            integer,          intent(in) :: nx, ny, nz
            double precision, intent(in) :: lx, ly, lz
            integer                      :: kx, ky, kz

            if (is_fft_initialised) then
                return
            endif

            is_fft_initialised = .true.

!             dz = dx(3)
!             dzi = dxi(3)
!             dz6  = f16 * dx(3)
!             dz2  = f12 * dx(3)
!             dz24 = f124 * dx(3)
!             dzisq = dxi(3) ** 2
!             hdzi = f12 * dxi(3)
!             nwx = nx / 2
!             nwy = ny / 2
!             nyp2 = ny + 2
!             nxp2 = nx + 2

            allocate(rkx(0:nx-1))
            allocate(hrkx(nx))
            allocate(rky(0:ny-1))
            allocate(hrky(ny))
            allocate(rkz(0:nz))
            allocate(rkzi(1:nz-1))
            allocate(xtrig(2 * nx))
            allocate(ytrig(2 * ny))
            allocate(ztrig(2 * nz))

            !----------------------------------------------------------------------
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, lx, ly, xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
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
            call init_deriv(nz, lz, rkz(1:nz))
            rkzi(1:nz-1) = one / rkz(1:nz-1)

        end subroutine init_fft

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        ! and returns the result as fs in spectral space (in x & y).
        ! Only FFTs over the x and y directions are performed.
        ! *** fp is destroyed upon exit ***
        subroutine fftxyp2s(fp, fs)
            double precision, intent(inout) :: fp(:, :, :)       !Physical
            double precision, intent(out)   :: fs(:, :, :)       !Spectral
            integer                         :: kx, iy, nzval, nxval, nyval

!             if (comm%size > 1) then
!                 call perform_fftxyp2s(fp, fs, xfactors, xtrig, yfactors, ytrig)
!             else

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
!             endif
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

            ! 1. Swap x and z from (z, x, y) to (x, z, y) pencil
            ! 2. Do x back-transform
            ! 3. Transform from (x, z, y) to (y, x, z) pencil
            ! 4. Do y back-transform
            ! 5. Transform from (y, x, z) to (z, y, x) pencil

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
        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine stzp2s(fs, ztrig, zfactors)
            double precision, intent(in) :: fs(box%hlo(3), box%hhi(3), &
                                               box%hlo(2), box%hhi(2), &
                                               box%hlo(1), box%hhi(1))
            double precision, intent(in) :: ztrig(:)
            integer,          intent(in) :: zfactors(5)
            integer                      :: kx, ky

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dst(1, nz, fs(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do
        end subroutine stzp2s

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ctzp2s(fs, ztrig, zfactors)
            double precision, intent(in) :: fs(box%hlo(3), box%hhi(3), &
                                               box%hlo(2), box%hhi(2), &
                                               box%hlo(1), box%hhi(1))
            double precision, intent(in) :: ztrig(:)
            integer,          intent(in) :: zfactors(5)
            integer                      :: kx, ky

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, fs(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

        end subroutine ctzp2s

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffx(fs, ds)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: ds(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision              :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            integer                       :: kx, dkx, kxc, nwx, nxp2
            double precision              :: si

            nwx = nx / 2
            nxp2 = nx + 1

            call reverse_x(fs, gs)

            !Carry out differentiation by wavenumber multiplication:
            if (0 == box%lo(1)) then
                ds(:, :, 0) = zero
            endif

            do kx = max(1, box%lo(1)), box%hi(1)
                dkx = min(2 * kx, 2 * (nx - kx))
                si = merge(1.0d0, -1.0d0, kx >= nwx + 1)
                ds(:, :, kx)  = si * hrkx(dkx) * gs(:, :, kx-1)
            enddo

            if (mod(nx, 2) .eq. 0) then
                kxc = nwx! + 1
                if (kxc >= box%lo(1) .and. kxc <= box%hi(1)) then
                    ds(:, :, kxc) = zero
                endif
            endif

        end subroutine diffx

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dy
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        subroutine diffy(fs, ds)
            double precision, intent(in)  :: fs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: ds(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision              :: gs(box%hlo(3):box%hhi(3), &
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            integer                       :: ky, dky, kyc, nwy, nyp2
            double precision              :: si

            nwy = ny / 2
            nyp2 = ny + 1

            call reverse_y(fs, gs)

            !Carry out differentiation by wavenumber multiplication:
            if (0 == box%lo(2)) then
                ds(:, 0, :) = zero
            endif

            do ky = max(1, box%lo(2)), box%hi(2)
                dky = min(2 * ky, 2 * (ny - ky))
                si = merge(1.0d0, -1.0d0, ky >= nwy + 1)
                ds(:, ky, :)  = si * hrky(dky) * gs(:, ky-1, :)
            enddo

            if (mod(ny, 2) .eq. 0) then
                kyc = nwy
                if (kyc >= box%lo(2) .and. kyc <= box%hi(2)) then
                    ds(:, kyc, :) = zero
                endif
            endif

        end subroutine diffy

end module sta3dfft
