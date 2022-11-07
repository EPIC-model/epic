! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_fft_3d_sphere
    use unit_test
    use constants, only : pi, twopi, f12, zero, five, four, three, two, one
    use inversion_utils, only : init_fft, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :), &
                                     fp2(:, :, :), &
                                     fp(:, :, :),  &
                                     fs(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z, rad


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
                ! offset by irrational numbers between 0 and 1
                rad = (x-sqrt(two)+one)**2 + (y-sqrt(three)+one)**2 + (z-sqrt(five)+two)**2
                if(rad<2.0) then 
                    fp1(k, j, i) = 1.+rad
                else
                    fp1(k, j, i) = 0.
                endif
                fp2(k, j, i) = fp1(k, j, i)
            enddo
        enddo
    enddo

    fs = zero

    ! we need to copy since fftxyp2s overwrites *fp*.
    fp = fp1

    ! forward FFT
    call fftxyp2s(fp, fs)

    ! inverse FFT
    call fftxys2p(fs, fp2)

    ! final check
    error = maxval(dabs(fp1 - fp2))

    call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-14))

end program test_fft_3d_sphere
