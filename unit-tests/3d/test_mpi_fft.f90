! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_fft_3d
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use inversion_utils, only : init_fft, xfactors, xtrig, yfactors, ytrig !, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use pencil_fft
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :), &
                                     fp2(:, :, :), &
                                     fs(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z
    logical                       :: passed

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(fp1(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fp2(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    call init_fft

    call initialise_pencil_fft(nx, ny, nz)

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

    ! we need to copy since fftxyp2s overwrites *fp*.

    ! forward FFT
    call perform_fftxyp2s(fp1(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          fs(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          xfactors, xtrig, yfactors, ytrig)

    ! inverse FFT
    call perform_fftxys2p(fs(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          fp2(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          xfactors, xtrig, yfactors, ytrig)
!     call fftxys2p(fs, fp2)

    call finalise_pencil_fft

    error = maxval(dabs(fp1(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)) - &
                        fp2(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        if (.not. passed) then
            call print_result_logical('Test FFT 2D transform', passed)
        else
            call print_result_dp('Test FFT 2D transform', error, atol=dble(1.0e-14))
        endif
    endif

end program test_fft_3d
