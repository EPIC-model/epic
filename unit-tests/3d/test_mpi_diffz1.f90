! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz using the
!  function:
!               cos(k * x) * sin(l * y) * cos(m * z)
!  where k = 2pi/L_x, l = 2pi/L_y and m = pi/L_z and where x, y and z all start
!  at 0 (one could start at -pi for x and y just as well).
!  The subroutine diffz should return
!               - m * cos(k * x) * sin(l * y) * sin(m * z)
! =============================================================================
program test_mpi_diffz1
    use unit_test
    use constants, only : zero, one, two, pi, twopi
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_mod, only : diffz, init_inversion, field_decompose_physical, field_combine_physical
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: fp(:, :, :), dp(:, :, :), &
                                     fs(:, :, :), ds(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z, k, l, m, prefactor
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    ! doubling nx, ny, nz reduces the error by a factor 16
    nx = 64
    ny = 64
    nz = 64
    lower  = (/zero, zero, zero/)
    extent =  (/pi, pi, pi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters


    allocate(fp(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(dp(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(ref_sol(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    k = twopi / extent(1)
    l = twopi / extent(2)
    m =    pi / extent(3)

    prefactor = - one / (k ** 2 + l ** 2 + m ** 2)

    dp = zero
    ref_sol = zero

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + iz * dx(3)
                fp(iz, iy, ix) = dcos(k * x) * dsin(l * y) * dcos(m * z)
                ref_sol(iz, iy, ix) = -m * dcos(k * x) * dsin(l * y) * dsin(m * z)
            enddo
        enddo
    enddo

    call init_inversion

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)

    error = maxval(dabs(dp(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                 - ref_sol(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

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

    passed = (passed .and. (comm%err == 0) .and. (error < 0.2d0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI diffz1', passed)
    endif

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

end program test_mpi_diffz1
