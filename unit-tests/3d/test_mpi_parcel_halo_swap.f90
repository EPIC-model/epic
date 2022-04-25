! =============================================================================
!                         Test MPI parcel halo swap
!
!   This unit test checks parcel halo swap. Each MPI rank sends mpi_rank+1
!   parcels to each of its neighbours.
! =============================================================================
program test_mpi_parcel_halo_swap
    use constants, only : zero, one, f12
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    use fields, only : field_alloc
    use parcel_container
    use parcel_mpi
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, dx
    use mpi_collectives
    use timer
    implicit none

    logical :: passed = .true.
    integer :: n_total, n, j, n_verify

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call update_parameters

    ! calls mpi_layout_init internally
    call field_alloc

    n_parcels = 8 * (mpi_rank + 1)
    n_total = 4 * mpi_size * (mpi_size + 1)
    call parcel_alloc(2 * n_total)

    do n = 1, mpi_rank + 1
        ! place parcels in southwest halo
        parcels%position(1, n) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n) = f12

        ! place parcels in west halo
        j = mpi_rank + 1
        parcels%position(1, n + j) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%lo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in northwest halo
        j = 2 * (mpi_rank + 1)
        parcels%position(1, n + j) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in north halo
        j = 3 * (mpi_rank + 1)
        parcels%position(1, n + j) = (box%lo(1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in northeast halo
        j = 4 * (mpi_rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in east halo
        j = 5 * (mpi_rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%lo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in southeast halo
        j = 6 * (mpi_rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in south halo
        j = 7 * (mpi_rank + 1)
        parcels%position(1, n + j) = (box%lo(1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12
    enddo

    parcels%volume(1:n_parcels) = mpi_rank + 1
    parcels%B(:, 1:n_parcels) = mpi_rank + 1
    parcels%vorticity(:, 1:n_parcels) = mpi_rank + 1
    parcels%buoyancy(1:n_parcels) = mpi_rank + 1

    call parcel_halo_swap

    n_verify = n_parcels
    call mpi_blocking_reduce(n_verify, MPI_SUM)

    passed = (passed .and. (n_verify == n_total))

    print *, n_total, n_verify

!     if (mpi_rank == mpi_master) then
!         total_vol = dble(n_total) * parcels%volume(1)
!         passed = (passed .and. (dabs(parcel_stats(IDX_KE) - 0.375d0 * total_vol) == zero))
!         passed = (passed .and. (dabs(parcel_stats(IDX_N_SMALL) - n_total) == zero))
!         passed = (passed .and. (dabs(parcel_stats(IDX_AVG_LAM) - one) < 1.0e-13))
!         passed = (passed .and. (parcel_stats(IDX_STD_LAM) < 1.0e-7))
!         passed = (passed .and. (parcel_stats(IDX_STD_VOL) < 1.0e-15))
!         passed = (passed .and. (dabs(parcel_stats(IDX_SUM_VOL) - total_vol) < 1.0e-15))
!         passed = (passed .and. (dabs(parcel_stats(IDX_NTOT_PAR) - n_total) == zero))
!         passed = (passed .and. (dabs(parcel_stats(IDX_RMS_XI) - f12) < 1.0e-15))
!         passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ETA) - f12) < 1.0e-15))
!         passed = (passed .and. (dabs(parcel_stats(IDX_RMS_ZETA) - f12) < 1.0e-15))
!     endif


    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI parcel halo swap', passed)
    endif

end program test_mpi_parcel_halo_swap
