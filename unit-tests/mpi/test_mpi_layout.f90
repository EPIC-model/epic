! =============================================================================
!                       Test MPI initialisation
!
!                   This unit test checks the MPI layout.
! =============================================================================
program test_mpi_layout
    use unit_test
    use mpi_communicator
    use mpi_layout
    implicit none

    integer, parameter            :: nx = 32, ny = 32, nz = 32
    double precision, parameter   :: lower(3) = (/0.0d0, 0.0d0, 0.0d0/)
    double precision, parameter   :: extent(3) = (/1.0d0, 1.0d0, 1.0d0/)
    double precision, allocatable :: data(:, :, :)
    double precision              :: sendbuf, recvbuf
    logical                       :: passed = .false.

    call mpi_comm_initialise

!     if (comm%size == 1) then
!         print *, "MPI tests must be run with more than one process."
!         call mpi_comm_finalise
!         stop
!     endif

    passed = (comm%err == 0)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    allocate(data(box%hlo(3):box%hhi(3), box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))

    data(:, :, :) = 1.0d0

    sendbuf = sum(data(box%lo(3):box%hi(3), box%lo(2):box%hi(2), box%lo(1):box%hi(1)))
    recvbuf = 0.0d0

    call MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm%world, comm%err)

    passed = (passed .and. (comm%err == 0) .and. (dble((nz+1)*nx*ny) - recvbuf == 0.0d0))

    call mpi_comm_finalise

    passed = (passed .and. (comm%err == 0))

    if (comm%rank == 0) then
        call print_result_logical('Test MPI layout', passed)
    endif

    deallocate(data)

end program test_mpi_layout
