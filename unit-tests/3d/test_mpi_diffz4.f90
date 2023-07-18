! =============================================================================
!                       Test subroutine diffz
!
!  This unit test checks the subroutine diffz and central_diffz using the
!  function:
!               z^2
!  in a domain of width 2 * pi in x and y, and of height pi (-pi/2 < z < pi/2)
!  The subroutine diffz and central_diffz should return
!               2z
! =============================================================================
program test_mpi_diffz4
    use unit_test
    use constants, only : zero, one, two, pi, twopi, f12
    use parameters, only : lower, update_parameters, dx, nx, ny, nz, extent
    use inversion_mod, only : init_inversion, central_diffz, diffz &
                            , field_combine_physical, field_decompose_physical, fftxyp2s, fftxys2p
    use mpi_communicator
    use mpi_layout
    implicit none

    double precision              :: error
    double precision, allocatable :: fs(:, :, :), fp(:, :, :), &
                                     ds(:, :, :), dp(:, :, :), &
                                     ref_sol(:, :, :)
    integer                       :: iz
    double precision              :: z
    logical                       :: passed = .false.

    call mpi_comm_initialise

    passed = (world%err == 0)

    nx = 128
    ny = 128
    nz = 64
    lower  = (/zero, zero, -f12 * pi/)
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

    do iz = 0, nz
        z = lower(3) + dble(iz) * dx(3)
        fp(iz, :, :) = z ** 2
        ref_sol(iz, :, :) = two * z
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

    call mpi_comm_finalise

    passed = (passed .and. (world%err == 0) .and. (error < 0.05d0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI diffz4', passed)
    endif

    deallocate(fp)
    deallocate(dp)
    deallocate(fs)
    deallocate(ds)
    deallocate(ref_sol)

end program test_mpi_diffz4
