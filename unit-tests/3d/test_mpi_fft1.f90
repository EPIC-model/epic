! =============================================================================
!                       Test FFT module
!
!       This unit test checks the FFT module with reference solutions.
! =============================================================================
program test_mpi_fft1
    use unit_test
    use constants, only : pi, twopi, f12, zero, five, four, three, two, one
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
    double precision              :: x, y, z, rad
    logical                       :: passed = .false.

    call mpi_env_initialise

    passed = (world%err == 0)

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

    ! forward FFT
    call fftxyp2s(fp1, fs)

    ! inverse FFT
    call fftxys2p(fs, fp2)

    call finalise_fft

    ! final check
    error = maxval(abs(fp1(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                      - fp2(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, world%root, world%comm, world%err)
    else
        call MPI_Reduce(error, error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, world%root, world%comm, world%err)
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0) .and. (error < dble(1.0e-14)))

    if (world%rank == world%root) then
        call print_result_logical('Test FFT 2D transform', passed)
    endif

    deallocate(fp1)
    deallocate(fp2)
    deallocate(fs)

end program test_mpi_fft1
