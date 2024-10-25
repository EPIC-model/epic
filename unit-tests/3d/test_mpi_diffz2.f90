! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz and central_diffz_semi_spectral
!  using the function:
!               cos(x) cos(y) cos(z)
!  in a domain of width 2 * pi in x and y, and of height pi (0 < z < pi)
!  The subroutine diffz and central_diffz_semi_spectral should return
!               - cos(x) cos(y) sin(z)
! =============================================================================
program test_mpi_diffz2
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_mod, only : central_diffz_semi_spectral, diffz, init_inversion, fftxyp2s, fftxys2p &
                            , field_decompose_physical, field_combine_physical
    use mpi_environment
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), fp(:, :, :), &
                                     ds(:, :, :), dp(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: ix, iy, iz
    double precision              :: x, y, z
    logical                       :: passed = .false.

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 128
    ny = 128
    nz = 64
    lower  = (/zero, zero, zero/)
    extent =  (/twopi, twopi, pi/)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

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
                fp(iz, iy, ix) = cos(x) * cos(y) * cos(z)
                ref_sol(iz, iy, ix) = -cos(x) * cos(y) * sin(z)
            enddo
        enddo
    enddo

    call init_inversion

    dp = fp
    call fftxyp2s(dp, fs)
    call central_diffz_semi_spectral(fs, ds)
    call fftxys2p(ds, dp)

    error = maxval(abs(dp(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                 - ref_sol(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

    call field_decompose_physical(fp, fs)
    call diffz(fs, ds)
    call field_combine_physical(ds, dp)


    error = max(error, maxval(abs(dp(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) &
                                 - ref_sol(:, box%lo(2):box%hi(2), box%lo(1):box%hi(1)))))

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

    passed = (passed .and. (world%err == 0) .and. (error < 4.0e-2))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI diffz2', passed)
    endif

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

    contains

end program test_mpi_diffz2
