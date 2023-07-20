! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_mpi_fft_3d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : initialise_fft, finalise_fft, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_environment
    use mpi_layout
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :), &
                                     fp2(:, :, :), &
                                     fs(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z
    logical                       :: passed

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    allocate(fp1(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fp2(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call initialise_fft(extent)

    fp1(:, :, :) = zero
    fp2(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                fp1(k, j, i) = dcos(four * x) + dsin(y) + dcos(z)
!                 fp2(k, j, i) = fp1(k, j, i)
            enddo
        enddo
    enddo

    fs = zero

    ! forward FFT
    call fftxyp2s(fp1, fs)

    ! inverse FFT
    call fftxys2p(fs, fp2)

    call finalise_fft

    error = maxval(dabs(fp1(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)) - &
                        fp2(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        if (.not. passed) then
            call print_result_logical('Test FFT 2D transform', passed)
        else
            call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-14))
        endif
    endif

end program test_mpi_fft_3d
