! =============================================================================
!                       Test horizontal spectral differentiation
!
!       This unit test checks the diffx and diffy subroutines.
! =============================================================================
program test_mpi_diffx
    use unit_test
    use constants, only : pi, twopi, f12, zero, four, two
    use parameters, only : update_parameters, dx, nx, ny, nz, lower, extent, upper
    use deriv1d
    use stafft
    use mpi_communicator
    use mpi_layout
    use fft_utils, only : reorder
    implicit none

    double precision              :: error = zero
    double precision, allocatable :: fp(:, :, :), &
                                     gp(:, :, :), &
                                     hp(:, :, :)
    integer                       :: i, j, k
    double precision              :: x, y, z, xr
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 8
    ny = 8
    nz = 8
    lower = (/-pi, -pi, -pi/)
!     extent = (/nx, ny, nz/)
    extent = (/twopi, twopi, twopi/)

    call update_parameters

    call mpi_layout_init(nx, ny, nz)

    allocate(fp(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(gp(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(hp(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

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
!                 fp(k, j, i) =          i + nx * j + nx * ny * k!dcos(four * x) + dsin(y) + dsin(z) * z
!                 gp(k, j, i) = nx - 1 - i + nx * j + nx * ny * k !fp(k, j, i)

!                 fp(k, j, i) = x + y + z
!                 gp(k, j, i) = xr + y + z
!                 print *, gp(k, j, i), fp(k, j, i)
                fp(k, j, i) = four * x + dsin(y) + dsin(z) * z
                gp(k, j, i) = four * xr + dsin(y) + dsin(z) * z
!                 print *, fp(k, j, i), gp(k, j, i)
            enddo
        enddo
    enddo

    call reorder(fp, hp)

    print *, "Reordering:", maxval(gp(box%lo(3):box%hi(3), &
                                      box%lo(2):box%hi(2), &
                                      box%lo(1):box%hi(1)) &
                                 - hp(box%lo(3):box%hi(3), &
                                      box%lo(2):box%hi(2), &
                                      box%lo(1):box%hi(1)))

!     print *, comm%rank, "gp:", gp(box%lo(3), box%lo(2), box%lo(1):box%hi(1))
!     print *, comm%rank, "hp:", hp(box%lo(3), box%lo(2), box%lo(1):box%hi(1))

   !  check result test field
    do i = box%lo(1), box%hi(1)
        do j = box%lo(2), box%hi(2)
            do k = box%lo(3), box%hi(3)
                print *, comm%rank, i, j, k, gp(k, j, i), hp(k, j, i), (gp(k, j, i) == hp(k, j, i))
                passed = (passed .and. (gp(k, j, i) - hp(k, j, i) < 1.0e-14))
            enddo
        enddo
    enddo


    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    print *, "passed", passed

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
