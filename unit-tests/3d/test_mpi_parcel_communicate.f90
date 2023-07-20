! =============================================================================
!                         Test MPI parcel halo swap
!
!   This unit test checks parcel halo swap. Each MPI rank sends cart%rank+1
!   parcels to each of its neighbours.
! =============================================================================
program test_mpi_parcel_communicate
    use constants, only : zero, one, f12
    use unit_test
    use mpi_environment
    use mpi_layout
    use fields, only : field_alloc
    use parcel_container
    use parcel_mpi
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, dx
    use mpi_collectives
    use mpi_timer
    implicit none

    logical :: passed = .true.
    integer :: n_total, n, j, n_total_verify
    integer :: n_local_verify, n_expected

    call mpi_env_initialise

    passed = (world%err == 0)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    ! calls mpi_layout_init internally
    call field_alloc

    n_parcels = 8 * (cart%rank + 1)
    n_total = 4 * cart%size * (cart%size + 1)
    call parcel_alloc(2 * n_total)

    n_total_parcels = n_total

    do n = 1, cart%rank + 1
        ! place parcels in southwest halo
        parcels%position(1, n) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n) = f12

        ! place parcels in west halo
        j = cart%rank + 1
        parcels%position(1, n + j) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%lo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in northwest halo
        j = 2 * (cart%rank + 1)
        parcels%position(1, n + j) = (box%hlo(1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in north halo
        j = 3 * (cart%rank + 1)
        parcels%position(1, n + j) = (box%lo(1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in northeast halo
        j = 4 * (cart%rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = ((box%hhi(2) - 1) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in east halo
        j = 5 * (cart%rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%lo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in southeast halo
        j = 6 * (cart%rank + 1)
        parcels%position(1, n + j) = ((box%hhi(1) - 1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12

        ! place parcels in south halo
        j = 7 * (cart%rank + 1)
        parcels%position(1, n + j) = (box%lo(1) + f12) * dx(1)
        parcels%position(2, n + j) = (box%hlo(2) + f12) * dx(2)
        parcels%position(3, n + j) = f12
    enddo

    parcels%volume(1:n_parcels) = cart%rank + 1
    parcels%B(:, 1:n_parcels) = cart%rank + 1
    parcels%vorticity(:, 1:n_parcels) = cart%rank + 1
    parcels%buoyancy(1:n_parcels) = cart%rank + 1

    call parcel_communicate

    n_total_verify = n_parcels
    call mpi_blocking_reduce(n_total_verify, MPI_SUM)

    ! this needs to be checked by the MPI root only!
    if (world%rank == world%root) then
        passed = (passed .and. (n_total_verify == n_total))
    endif

    n_local_verify = (neighbours(MPI_WEST)%rank+1)      &
                   + (neighbours(MPI_EAST)%rank+1)      &
                   + (neighbours(MPI_NORTH)%rank+1)     &
                   + (neighbours(MPI_SOUTH)%rank+1)     &
                   + (neighbours(MPI_NORTHEAST)%rank+1) &
                   + (neighbours(MPI_NORTHWEST)%rank+1) &
                   + (neighbours(MPI_SOUTHEAST)%rank+1) &
                   + (neighbours(MPI_SOUTHWEST)%rank+1)


    passed = (passed .and. (n_parcels == n_local_verify))

    if (passed) then
        n_total = int(sum(parcels%volume(1:n_parcels)))

        n_expected = (neighbours(MPI_WEST)%rank+1) ** 2      &
                   + (neighbours(MPI_EAST)%rank+1) ** 2      &
                   + (neighbours(MPI_NORTH)%rank+1) ** 2     &
                   + (neighbours(MPI_SOUTH)%rank+1) ** 2     &
                   + (neighbours(MPI_NORTHEAST)%rank+1) ** 2 &
                   + (neighbours(MPI_NORTHWEST)%rank+1) ** 2 &
                   + (neighbours(MPI_SOUTHEAST)%rank+1) ** 2 &
                   + (neighbours(MPI_SOUTHWEST)%rank+1) ** 2

        passed = (passed .and. (n_expected == n_total))
    endif

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel halo swap', passed)
    endif

    call mpi_env_finalise

end program test_mpi_parcel_communicate
