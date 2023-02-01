! =============================================================================
!                       Test horizontal spectral differentiation
!
!       This unit test checks the diffx and diffy subroutines.
! =============================================================================
program test_mpi_diffy
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use inversion_utils, only : init_fft, xfactors, xtrig, yfactors, ytrig
!     use inversion_utils, only : init_fft, fftxyp2s, fftxys2p
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent
    use deriv1d
    use stafft
    use fft_utils, only : diffy
    use pencil_fft
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

    nx = 128
    ny = 64
    nz = 64
    lower = (/-pi, -pi, -pi/)
    extent = (/twopi, twopi, twopi/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(fp(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(ds(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    call init_fft

    call initialise_pencil_fft(nx, ny, nz)

    fp(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                fp(k, j, i) = dcos(four * y)
            enddo
        enddo
    enddo

    fs = zero

    ! forward FFT
    call perform_fftxyp2s(fp(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          fs(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          xfactors, xtrig, yfactors, ytrig)

    fp = zero

    call diffy(fs, ds)

    ! inverse FFT
    call perform_fftxys2p(ds(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          fp(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)), &
                          xfactors, xtrig, yfactors, ytrig)

    ! check result test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                passed = (passed .and. (fp(k, j, i) - (-four * dsin(four * y)) < 1.0e-12))
            enddo
        enddo
    enddo

    do i = box%lo(2), box%hi(2)
        y = lower(2) + dble(i) * dx(2)
        print *, comm%rank, y, fp(box%lo(3), i, box%lo(1)), -four * dsin(four * y)
    enddo

    call finalise_pencil_fft

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
            call print_result_logical('Test diffy', passed)
        else
            call print_result_dp('Test diffy', error, atol=dble(1.0e-14))
        endif
    endif

end program test_mpi_diffy
