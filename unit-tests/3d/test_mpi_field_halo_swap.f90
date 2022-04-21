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

    integer, parameter            :: nx = 4, ny = 4, nz = 1
    double precision, allocatable :: values(:, :, :)
    logical                       :: passed = .true.
    double precision              :: diff

    call mpi_comm_initialise

    passed = (mpi_err == 0)

    call mpi_layout_init(nx, ny, nz)

    allocate(values(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    ! set all values to 0.5
    values(:, :, :) = f12


    ! set interior values to 1, excluding overlap with halo region of neighbours
    values(box%lo(3):box%hi(3), box%lo(2)+2:box%hi(2)-1, box%lo(1)+2:box%hi(1)-1) = one

!     ! set halo corners to 1/4
!     values(box%lo(3):box%hi(3), box%hlo(2), box%hlo(1)) = f14
!     values(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hlo(1)) = f14
!     values(box%lo(3):box%hi(3), box%hhi(2)-1:box%hhi(2), box%hhi(1)-1:box%hhi(1)) = f14
!     values(box%lo(3):box%hi(3), box%hlo(2), box%hhi(1)-1:box%hhi(1)) = f14

    ! set interior corners to 1/4
!     values(box%lo(3):box%hi(3), box%hi(2), box%hi(1)) = f14
!     values(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%hi(1)) = f14
!     values(box%lo(3):box%hi(3), box%lo(2):box%lo(2)+1, box%lo(1):box%lo(1)+1) = f14
!     values(box%lo(3):box%hi(3), box%hi(2), box%lo(1):box%lo(1)+1) = f14

    call field_halo_swap(values)

    if (mpi_rank == 0) then
        print *, values(0, :, :)
!         print *, values(2, :)
! !         print *, values(, :)
!         print *, values(1, :)
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
