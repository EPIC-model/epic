! =============================================================================
!                       Test horizontal spectral differentiation
!
!       This unit test checks the diffx and diffy subroutines.
! =============================================================================
program test_mpi_diffx
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use sta3dfft, only : initialise_fft, finalise_fft, diffx, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp(:, :, :), &
                                     fs(:, :, :), &
                                     ds(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 64
    ny = 128
    nz = 64
    lower = (/-pi, -pi, -pi/)
    extent = (/twopi, twopi, twopi/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(fp(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(box%lo(3):box%hi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(ds(box%lo(3):box%hi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    call initialise_fft(extent)

    fp(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                fp(k, j, i) = dcos(four * x)
            enddo
        enddo
    enddo

    fs = zero

    ! forward FFT
    call fftxyp2s(fp, fs)

    fp = zero

    call diffx(fs, ds)

    ! inverse FFT
    call fftxys2p(ds, fp)

    ! check result test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                passed = (passed .and. (fp(k, j, i) - (-four * dsin(four * x)) < 1.0e-12))
            enddo
        enddo
    enddo

    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        print *, comm%rank, x, fp(box%lo(3), box%lo(2), i), -four * dsin(four * x)
    enddo

    call finalise_fft

    print *, "Layout:", layout%l_parallel

!     error = maxval(dabs(fp - dr))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

!     print *, error

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        if (.not. passed) then
            call print_result_logical('Test diffx', passed)
        else
            call print_result_dp('Test diffx', error, atol=dble(1.0e-14))
        endif
    endif

end program test_mpi_diffx
