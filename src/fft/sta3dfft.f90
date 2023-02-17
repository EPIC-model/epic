module sta3dfft
    use mpi_layout
    use mpi_reverse, only : reverse_x, reverse_y
    use stafft, only : dct, dst
    use constants, only : zero, one
    use stafft
    use sta2dfft
    use deriv1d, only : init_deriv
    use fft_pencil
    implicit none

    private

    ! Wavenumbers:
    double precision, protected, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:), rkzi(:)

    !Quantities needed in FFTs:
    double precision, protected, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer,          protected              :: xfactors(5), yfactors(5), zfactors(5)

    integer :: nwx, nwy

    integer :: nx, ny, nz

    logical :: is_fft_initialised = .false.

    public :: initialise_fft &
            , finalise_fft   &
            , diffx          &
            , diffy          &
            , fftxyp2s       &
            , fftxys2p       &
            , fftsine        &
            , fftcosine      &
            , rkx            &
            , rky            &
            , rkz            &
            , rkzi           &
            , zfactors       &
            , ztrig          &
            , xfactors       &
            , xtrig          &
            , yfactors       &
            , ytrig

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine initialise_fft(extent)
            double precision, intent(in) :: extent(3)
            integer                      :: kx, ky!, kz

            if (is_fft_initialised) then
                return
            endif

            is_fft_initialised = .true.

            nx = box%global_size(1)
            ny = box%global_size(2)
            nz = box%global_size(3)

            call initialise_pencil_fft(nx, ny, nz)

            nwx = nx / 2
            nwy = ny / 2

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
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(0) = zero
            do kx = 1, nwx-1
                rkx(kx)    = hrkx(2 * kx)
                rkx(nx-kx) = hrkx(2 * kx)
            enddo
            rkx(nwx) = hrkx(nx)

            !Define y wavenumbers:
            rky(0) = zero
            do ky = 1, nwy-1
                rky(ky)    = hrky(2 * ky)
                rky(ny-ky) = hrky(2 * ky)
            enddo
            rky(nwy) = hrky(ny)

            !Define z wavenumbers:
            rkz(0) = zero
            call init_deriv(nz, extent(3), rkz(1:nz))
            rkzi(1:nz-1) = one / rkz(1:nz-1)

        end subroutine initialise_fft

        subroutine finalise_fft
            if (allocated(rkx)) then
                deallocate(rkx)
                deallocate(hrkx)
                deallocate(rky)
                deallocate(hrky)
                deallocate(rkz)
                deallocate(rkzi)
                deallocate(xtrig)
                deallocate(ytrig)
                deallocate(ztrig)
            endif

            call finalise_pencil_fft
        end subroutine finalise_fft

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes a 2D FFT (in x & y) of a 3D array fp in physical space
        ! and returns the result as fs in spectral space (in x & y).
        ! Only FFTs over the x and y directions are performed.
        ! *** fp is destroyed upon exit ***
        subroutine fftxyp2s(fp, fs)
            double precision, intent(in)  :: fp(box%hlo(3):box%hhi(3), & !Physical
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            double precision, intent(out) :: fs(box%lo(3):box%hi(3),   & !Spectral
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            integer                       :: i, j

            ! 1. Transform from (z, y, x) to (y, x, z) pencil
            ! 2. Do y transform
            ! 3. Transform from (y, x, z) to (x, z, y) pencil
            ! 4. Do x transform
            ! 5. Transform from (x, z, y) to (y, x, z) pencil
            ! 6. Transform from (y, x, z) to (z, y, x) pencil

            call transpose_to_pencil(y_from_z_transposition,  &
                                     (/1, 2, 3/),             &
                                     dim_y_comm,              &
                                     FORWARD,                 &
                                     fp(box%lo(3):box%hi(3),  &
                                        box%lo(2):box%hi(2),  &
                                        box%lo(1):box%hi(1)), &
                                     fft_in_y_buffer)

            do i = 1, size(fft_in_y_buffer, 2)
                do j = 1, size(fft_in_y_buffer, 3)
                    call forfft(1, size(fft_in_y_buffer, 1), fft_in_y_buffer(:, i, j), ytrig, yfactors)
                enddo
            enddo

            call transpose_to_pencil(x_from_y_transposition, &
                                     (/2, 3, 1/),            &
                                     dim_x_comm,             &
                                     FORWARD,                &
                                     fft_in_y_buffer,        &
                                     fft_in_x_buffer)

            do i = 1, size(fft_in_x_buffer, 2)
                do j = 1, size(fft_in_x_buffer, 3)
                    call forfft(1, size(fft_in_x_buffer, 1), fft_in_x_buffer(:, i, j), xtrig, xfactors)
                enddo
            enddo

            call transpose_to_pencil(y_from_x_transposition,    &
                                     (/3, 1, 2/),               &
                                     dim_x_comm,                &
                                     BACKWARD,                  &
                                     fft_in_x_buffer,           &
                                     fft_in_y_buffer)

            call transpose_to_pencil(z_from_y_transposition,  &
                                     (/2, 3, 1/),             &
                                     dim_y_comm,              &
                                     BACKWARD,                &
                                     fft_in_y_buffer,         &
                                     fs)

        end subroutine fftxyp2s

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Computes an *inverse* 2D FFT (in x & y) of a 3D array fs in spectral
        ! space and returns the result as fp in physical space (in x & y).
        ! Only inverse FFTs over the x and y directions are performed.
        ! *** fs is destroyed upon exit ***
        subroutine fftxys2p(fs, fp)
            double precision, intent(in)  :: fs(box%lo(3):box%hi(3),   & !Spectral
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision, intent(out) :: fp(box%hlo(3):box%hhi(3), & !Physical
                                                box%hlo(2):box%hhi(2), &
                                                box%hlo(1):box%hhi(1))
            integer                       :: i, j

            ! 1. Transform to (z, y, x) to (y, x, z) pencil
            ! 2. Do y back-transform
            ! 3. Transform from (y, x, z) to (x, z, y) pencil
            ! 4. Do x back-transform
            ! 5. Transform from (x, z, y) to (y, x, z) pencil
            ! 6. Transform from (y, x, z) to (z, y, x) pencil

            call transpose_to_pencil(y_from_z_transposition, &
                                     (/1, 2, 3/),            &
                                     dim_y_comm,             &
                                     FORWARD,                &
                                     fs,                     &
                                     fft_in_y_buffer)

            call transpose_to_pencil(x_from_y_transposition, &
                                     (/2, 3, 1/),            &
                                     dim_x_comm,             &
                                     FORWARD,                &
                                     fft_in_y_buffer,        &
                                     fft_in_x_buffer)

            do i = 1, size(fft_in_x_buffer, 2)
                do j = 1, size(fft_in_x_buffer, 3)
                    call revfft(1, size(fft_in_x_buffer, 1), fft_in_x_buffer(:, i, j), xtrig, xfactors)
                enddo
            enddo

            call transpose_to_pencil(y_from_x_transposition, &
                                     (/3, 1, 2/),            &
                                     dim_x_comm,             &
                                     BACKWARD,               &
                                     fft_in_x_buffer,        &
                                     fft_in_y_buffer)

            do i = 1, size(fft_in_y_buffer, 2)
                do j = 1, size(fft_in_y_buffer, 3)
                    call revfft(1, size(fft_in_y_buffer, 1), fft_in_y_buffer(:, i, j), ytrig, yfactors)
                enddo
            enddo

            call transpose_to_pencil(z_from_y_transposition,  &
                                     (/2, 3, 1/),             &
                                     dim_y_comm,              &
                                     BACKWARD,                &
                                     fft_in_y_buffer,         &
                                     fp(box%lo(3):box%hi(3),  &
                                        box%lo(2):box%hi(2),  &
                                        box%lo(1):box%hi(1)))

        end subroutine fftxys2p

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine fftsine(fs)
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), &  ! 0:nz
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                         :: kx, ky

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
!                     fs(1:nz, ky, kx) = zero !FIXME ask David
                    call dst(1, nz, fs(1:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do
        end subroutine fftsine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine fftcosine(fs)
            double precision, intent(inout) :: fs(box%lo(3):box%hi(3), & ! 0:nz
                                                  box%lo(2):box%hi(2), &
                                                  box%lo(1):box%hi(1))
            integer                         :: kx, ky

            !$omp parallel do collapse(2) private(kx, ky)
            do kx = box%lo(1), box%hi(1)
                do ky = box%lo(2), box%hi(2)
                    call dct(1, nz, fs(0:nz, ky, kx), ztrig, zfactors)
                enddo
            enddo
            !$omp end parallel do

        end subroutine fftcosine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dx
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        ! Note: gs must have halo in x due to the reordering algorithm.
        subroutine diffx(fs, ds)
            double precision, intent(in)  :: fs(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision              :: gs(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%lo(2):box%hi(2),   &
                                                box%hlo(1):box%hhi(1))
            integer                       :: kx, dkx
            double precision              :: si


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
                if (nwx >= box%lo(1) .and. nwx <= box%hi(1)) then
                    ds(:, :, nwx) = zero
                endif
            endif

        end subroutine diffx

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Given fs in spectral space (at least in x & y), this returns dfs/dy
        ! (partial derivative).  The result is returned in ds, again
        ! spectral.  Uses exact form of the derivative in spectral space.
        ! Note: gs must have halo in y due to the reordering algorithm.
        subroutine diffy(fs, ds)
            double precision, intent(in)  :: fs(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision, intent(out) :: ds(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%lo(2):box%hi(2),   &
                                                box%lo(1):box%hi(1))
            double precision              :: gs(box%lo(3):box%hi(3),   & ! 0:nz
                                                box%hlo(2):box%hhi(2), &
                                                box%lo(1):box%hi(1))
            integer                       :: ky, dky
            double precision              :: si

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
                if (nwy >= box%lo(2) .and. nwy <= box%hi(2)) then
                    ds(:, nwy, :) = zero
                endif
            endif

        end subroutine diffy

end module sta3dfft
