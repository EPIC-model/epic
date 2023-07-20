! =============================================================================
!                       Test MPI reverse_y subroutine
!
!       This unit test checks the reverse algorithm used in diffy.
! =============================================================================
program test_mpi_reverse_y
    use unit_test
    use constants, only : zero, one
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent, upper
    use deriv1d
    use stafft
    use mpi_environment
    use mpi_layout
    use mpi_reverse, only : reverse_y               &
                          , intialise_mpi_reverse   &
                          , finalise_mpi_reverse
    implicit none

    double precision, allocatable :: fp(:, :, :), &
                                     gp(:, :, :), &
                                     hp(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z, xr
    logical                       :: passed = .false.

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 8
    ny = 8
    nz = 8
    lower = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call intialise_mpi_reverse

    allocate(fp(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(gp(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%hlo(1):box%hhi(1)))
    allocate(hp(box%lo(3):box%hi(3), box%hlo(2):box%hhi(2), box%lo(1):box%hi(1)))

    fp(:, :, :) = zero
    gp(:, :, :) = zero
    hp(:, :, :) = zero

    ! setup test field
    do i = box%lo(1), box%hi(1)
        x = lower(1) + dble(i) * dx(1)
        xr = upper(1) - dx(1) - dble(i) * dx(1)
        do j = box%lo(2), box%hi(2)
            y = lower(2) + dble(j) * dx(2)
            do k = box%lo(3), box%hi(3)
                z = lower(3) + dble(k) * dx(3)
                fp(k, j, i) = i + nx *            j + nx * ny * k
                gp(k, j, i) = i + nx * (ny - 1 - j) + nx * ny * k
            enddo
        enddo
    enddo

    call reverse_y(fp, hp)

   !  check result test field
    do i = box%lo(1), box%hi(1)
        do j = box%lo(2), box%hi(2)
            do k = box%lo(3), box%hi(3)
                passed = (passed .and. (gp(k, j, i) - hp(k, j, i) < 1.0e-14))
            enddo
        enddo
    enddo

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    call finalise_mpi_reverse

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
            call print_result_logical('Test MPI reverse_y', passed)
    endif

end program test_mpi_reverse_y
