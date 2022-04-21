! =============================================================================
!                             Test MPI field halo
!
!                   This unit test checks filling the field halo.
! =============================================================================
program test_mpi_layout
    use constants, only : zero
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    implicit none

    integer, parameter            :: nx = 8, ny = 10, nz = 4
    double precision, allocatable :: values(:, :, :)
    logical                       :: passed = .true.
    double precision              :: diff

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    values(:, :, :) = dble(mpi_rank + 1)

    call field_halo_fill(values)


    !
    ! check results
    !

    ! check west halo
    diff = sum(abs(dble(neighbour%west + 1) - values(box%lo(3):box%hi(3), &
                                                     box%lo(2):box%hi(2), &
                                                     box%hlo(1))))

    ! check the east halo
    diff = diff + sum(abs(dble(neighbour%east + 1) - values(box%lo(3):box%hi(3),       &
                                                            box%lo(2):box%hi(2),       &
                                                            box%hhi(1)-1:box%hhi(1))))

    ! check the south halo
    diff = diff + sum(abs(dble(neighbour%south + 1) - values(box%lo(3):box%hi(3),      &
                                                             box%hlo(2),               &
                                                             box%lo(1):box%hi(1))))

    ! check the north halo
    diff = diff + sum(abs(dble(neighbour%north + 1) - values(box%lo(3):box%hi(3),      &
                                                             box%hhi(2)-1:box%hhi(2),  &
                                                             box%lo(1):box%hi(1))))

    ! check the southwest halo
    diff = diff + sum(abs(dble(neighbour%southwest + 1) - values(box%lo(3):box%hi(3),  &
                                                                 box%hlo(2),           &
                                                                 box%hlo(1))))

    ! check the northwest halo
    diff = diff + sum(abs(dble(neighbour%northwest + 1) - values(box%lo(3):box%hi(3),       &
                                                                 box%hhi(2)-1:box%hhi(2),   &
                                                                 box%hlo(1))))

    ! check the northeast halo
    diff = diff + sum(abs(dble(neighbour%northeast + 1) - values(box%lo(3):box%hi(3),       &
                                                                 box%hhi(2)-1:box%hhi(2),   &
                                                                 box%hhi(1)-1:box%hhi(1))))

    ! check the southeast halo
    diff = diff + sum(abs(dble(neighbour%southeast + 1) - values(box%lo(3):box%hi(3),       &
                                                                 box%hlo(2),                &
                                                                 box%hhi(1)-1:box%hhi(1))))

    passed = (passed .and. (diff == zero))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (mpi_err == 0))

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI field halo', passed)
    endif

    deallocate(values)

end program test_mpi_layout
