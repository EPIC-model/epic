! =============================================================================
!                             Test MPI field halo swap
!
!                   This unit test checks filling interior and halo.
! =============================================================================
program test_field_halo_swap
    use constants, only : zero, f12, f14, one
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    implicit none

    integer, parameter            :: nx = 6, ny = 6, nz = 0
    double precision, allocatable :: values(:, :, :)
    logical                       :: passed = .true.
    double precision              :: diff
    integer                       :: i, j

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    ! set all values to 1
    values(:, :, :) = one


    call field_halo_swap(values)

    if (mpi_rank == 0) then
        do i = box%hlo(1), box%hhi(1)
            do j = box%hlo(2), box%hhi(2)
                print *, j, i, values(0, j, i)
            enddo
        enddo
    endif

    call mpi_comm_finalise
    stop

    !
    ! check results
    !

    ! check west halo
    diff = sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    print *, diff

    ! check the east halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the south halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the north halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the southwest halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the northwest halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the northeast halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    ! check the southeast halo
    diff = diff + sum(abs(one - values(box%lo(3):box%hi(3), :, :)))

    passed = (passed .and. (diff == zero))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI field halo swap', passed)
    endif

    deallocate(values)

end program test_field_halo_swap
