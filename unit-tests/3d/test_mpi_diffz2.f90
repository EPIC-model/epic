! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz and central_diffz using the
!  function:
!               cos(x) cos(y) cos(z)
!  in a domain of width 2 * pi in x and y, and of height pi (0 < z < pi)
!  The subroutine diffz and central_diffz should return
!               - cos(x) cos(y) sin(z)
! =============================================================================
program test_mpi_diffz2
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_mod, only : central_diffz, diffz, init_inversion, fftxyp2s, fftxys2p &
                            , field_decompose_physical, field_combine_physical
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), fp(:, :, :), &
                                     ds(:, :, :), dp(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (comm%err == 0)

    nx = 128
    ny = 128
    nz = 64
    lower  = (/zero, zero, zero/)
    extent =  (/twopi, twopi, pi/)

    call update_parameters

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(fp(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(dp(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
    allocate(fs(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(ds(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    allocate(ref_sol(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    dp = zero
    ref_sol = zero

    do ix = box%lo(1), box%hi(1)
        x = lower(1) + ix * dx(1)
        do iy = box%lo(2), box%hi(2)
            y = lower(2) + iy * dx(2)
            do iz = 0, nz
                z = lower(3) + dble(iz) * dx(3)
                fp(iz, iy, ix) = dcos(x) * dcos(y) * dcos(z)
                ref_sol(iz, iy, ix) = -dcos(x) * dcos(y) * dsin(z)
            enddo
        enddo
    enddo

    call init_inversion

    dp = fp
    call fftxyp2s(dp, fs)
    call central_diffz(fs, ds)
    call fftxys2p(ds, dp)

    error = maxval(dabs(dp(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                 - ref_sol(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)


    error = max(error, maxval(dabs(dp(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                                 - ref_sol(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))))

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

    passed = (passed .and. (comm%err == 0) .and. (error < 4.0e-2))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI diffz2', passed)
    endif

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

    contains

end program test_mpi_diffz2
