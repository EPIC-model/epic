! =============================================================================
!                       Test MPI field interior accumulate
!
!   This unit test checks accumulating the field interior. Each grid point
!   (including halo) is assigned the value of 1. After
!   accumulating the interior values, each interior grid point overlapping
!   with halo grid points from other processes should be assigned a value of
!   n+1, where n is the number of processes that share the grid point.
! =============================================================================
program test_field_interior_accumulate
    use constants, only : zero, one, two, four
    use unit_test
    use mpi_communicator
    use mpi_layout
    use field_mpi
    implicit none

    integer, parameter            :: nx = 10, ny = 10, nz = 1
    double precision, parameter   :: lower(3) = (/zero, zero, zero/)
    double precision, parameter   :: extent(3) = (/one, one, one/)
    double precision, allocatable :: values(:, :, :)
    logical                       :: passed = .true.
    double precision              :: diff

    call mpi_comm_initialise

    if (comm%size > 4) then
        if (comm%rank == comm%master) then
            print *, "This unit test does not work with more than 4 MPI ranks."
        endif
        call mpi_comm_finalise
        stop
    endif

    passed = (comm%err == 0)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    values(:, :, :) = one

    call field_interior_accumulate_scalar(values, l_alloc=.true.)

    !
    ! check results
    !

    ! check non-overlapping interior grid points; they should still be 1
    diff = sum(abs(one - values(box%lo(3):box%hi(3), box%lo(2)+2:box%hi(2)-1, box%lo(1)+2:box%hi(1)-1)))

    ! check the southwest interior (i.e. corner); should have value 4
    diff = diff + sum(abs(four - values(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1)))

    ! check the northwest interior (i.e. corner); should have value 4
    diff = diff + sum(abs(four - values(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%lo(1)+1)))

    ! check the northeast interior (i.e. corner); should have value 4
    diff = diff + sum(abs(four - values(box%lo(3):box%hi(3), box%hi(2), box%hi(1))))

    ! check the southeast interior (i.e. corner); should have value 4
    diff = diff + sum(abs(four - values(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1))))

    ! check the west interior (excluding corners); should have values 2
    diff = diff + sum(abs(two - values(box%lo(3):box%hi(3), box%lo(2)+2:box%hi(2)-1, box%lo(1):box%lo(1)+1)))

    ! check the east interior (excluding corners); should have values 2
    diff = diff + sum(abs(two - values(box%lo(3):box%hi(3), box%lo(2)+2:box%hi(2)-1, box%hi(1))))

    ! check the south interior (excluding corners); should have values 2
    diff = diff + sum(abs(two - values(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1)+2:box%hi(1)-1)))

    ! check the north interior (excluding corners); should have values 2
    diff = diff + sum(abs(two - values(box%lo(3):box%hi(3), box%hi(2), box%lo(1)+2:box%hi(1)-1)))

    passed = (passed .and. (diff == zero))

    if (comm%rank == comm%master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, comm%master, comm%world, comm%err)
    endif

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == comm%master) then
        call print_result_logical('Test MPI field interior accumulate', passed)
    endif

    deallocate(values)

end program test_field_interior_accumulate
