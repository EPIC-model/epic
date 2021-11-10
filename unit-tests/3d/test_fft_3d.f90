! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_fft_3d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use inversion_utils_mod, only : init_fft, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :), &
                                     fp2(:, :, :), &
                                     fs(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z


    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call update_parameters

    allocate(fp1(0:nz, 0:ny-1, 0:nx-1))
    allocate(fp2(0:nz, 0:ny-1, 0:nx-1))
    allocate(fs(0:nz, 0:nx-1, 0:ny-1))

    call init_fft

    ! setup test field
    do i = 0, nx-1
        x = lower(1) + dble(i) * dx(1)
        do j = 0, ny-1
            y = lower(2) + dble(j) * dx(2)
            do k = 0, nz
                z = lower(3) + dble(k) * dx(3)
                fp1(k, j, i) = dcos(four * x) + dsin(y) + dcos(z)
                fp2(k, j, i) = fp1(k, j, i)
            enddo
        enddo
    enddo

    fs = zero

    ! forward FFT
    call fftxyp2s(fp1, fs)

    ! inverse FFT
    call fftxys2p(fs, fp2)

    ! final check
    error = maxval(dabs(fp1 - fp2))

    call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-14))

end program test_fft_3d
