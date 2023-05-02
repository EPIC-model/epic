! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_mpi_fft0
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : initialise_fft, finalise_fft, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp1(:, :, :), &
                                     fp2(:, :, :), &
                                     fs(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 32
    ny = 64
    nz = 128
    lower = (/-pi, f12 * pi, zero/)
    extent = (/twopi, two * twopi, four * twopi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    allocate(fp1(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fp2(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))

    call initialise_fft(extent)

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
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

    call finalise_fft

    ! final check
    error = maxval(dabs(fp1(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                      - fp2(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(error, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%master, comm%world, comm%err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0) .and. (error < dble(1.0e-14)))

    if (comm%rank == comm%master) then
        call print_result_logical('Test FFT 2D transform', passed)
    endif

    deallocate(fp1)
    deallocate(fp2)
    deallocate(fs)

end program test_mpi_fft0
