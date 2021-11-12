! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module for combined inner
!       dimensions.
! =============================================================================
program test_ft_2_3d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use inversion_utils, only : init_fft &
                              , xfactors &
                              , yfactors &
                              , xtrig    &
                              , ytrig
    use sta2dfft
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fpp(:, :, :), &
                                     fs1(:, :, :), &
                                     fs2(:, :, :), &
                                     fp1(:, :, :), &
                                     fp2(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z


    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call update_parameters

    allocate(fpp(-1:nz+1, 0:ny-1, 0:nx-1))
    allocate(fs1(0:nz, 0:nx-1, 0:ny-1))
    allocate(fs2(0:nz, 0:nx-1, 0:ny-1))
    allocate(fp1(-1:nz+1, 0:ny-1, 0:nx-1))
    allocate(fp2(-1:nz+1, 0:ny-1, 0:nx-1))

    call init_fft

    fpp = zero

    ! setup test field
    do i = 0, nx-1
        x = lower(1) + dble(i) * dx(1)
        do j = 0, ny-1
            y = lower(2) + dble(j) * dx(2)
            do k = 0, nz
                z = lower(3) + dble(k) * dx(3)
                fpp(k, j, i) = dcos(four * x) + dsin(y) + dcos(z)
            enddo
        enddo
    enddo

    fs1 = zero
    fs2 = zero

    fp1 = zero
    fp2 = zero

    ! forward FFT
    call combined_forward_fft(fpp(0:nz, :, :), fs1)

    call not_combined_forward_fft(fpp(0:nz, :, :), fs2)

    error = maxval(dabs(fs1 - fs2))

    ! backward FFT
    call combined_backward_fft(fs1, fp1(0:nz, :, :))

    call not_combined_backward_fft(fs1, fp2(0:nz, :, :))

    error = maxval(dabs(fp1 - fp2))

    call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-15))

    contains
        subroutine combined_forward_fft(fp, fs)
            double precision, intent(in)  :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out) :: fs(0:nz, nx, ny)  !Spectral
            double precision              :: fp_copy(0:nz, ny, nx)
            integer:: kx, iy

            fp_copy = fp

             !Carry out a full x transform first:
            call forfft((nz+1) * ny, nx, fp_copy, xtrig, xfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fs(:, kx, iy) = fp_copy(:, iy, kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            call forfft((nz+1)*nx, ny, fs, ytrig, yfactors)
        end subroutine combined_forward_fft

        subroutine not_combined_forward_fft(fp, fs)
            double precision, intent(in)  :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out) :: fs(0:nz, nx, ny)  !Spectral
            double precision              :: fp_copy(0:nz, ny, nx)
            integer                       :: kx, iy, iz

            fp_copy = fp

            !Carry out a full x transform first:
            do iz = 0, nz
                call forfft(ny, nx, fp_copy(iz,  :,  :), xtrig, xfactors)
            enddo

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fs(:, kx, iy) = fp_copy(:, iy, kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            do iz = 0, nz
                call forfft(nx, ny, fs(iz,  :,  :), ytrig, yfactors)
            enddo
        end subroutine not_combined_forward_fft


        subroutine combined_backward_fft(fs, fp)
            double precision, intent(in)  :: fs(0:nz, nx, ny)  !Spectral
            double precision, intent(out) :: fp(0:nz, ny, nx)  !Physical
            double precision              :: fs_copy(0:nz, nx, ny)
            integer                       :: kx, iy

            fs_copy = fs

            !Carry out a full inverse y transform first:
            call revfft((nz+1) * nx, ny, fs_copy, ytrig, yfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fp(:, iy, kx) = fs_copy(:, kx, iy)
                enddo
            enddo

            !Carry out a full inverse x transform:
            call revfft((nz+1) * ny, nx, fp, xtrig, xfactors)
        end subroutine combined_backward_fft


        subroutine not_combined_backward_fft(fs, fp)
            double precision, intent(in)  :: fs(0:nz, nx, ny)  !Spectral
            double precision, intent(out) :: fp(0:nz, ny, nx)  !Physical
            double precision              :: fs_copy(0:nz, nx, ny)
            integer                       :: kx, iy, iz

            fs_copy = fs

            !Carry out a full inverse y transform first:
            do iz = 0, nz
                call revfft(nx, ny, fs_copy(iz, :, :), ytrig, yfactors)
            enddo

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fp(:, iy, kx) = fs_copy(:, kx, iy)
                enddo
            enddo

            !Carry out a full inverse x transform:
            do iz = 0, nz
                call revfft(ny, nx, fp(iz,  :,  :), xtrig, xfactors)
            enddo
        end subroutine not_combined_backward_fft

end program test_ft_2_3d
