! =============================================================================
!                       Test MPI initialisation
!
!               This unit test checks the MPI init and finalize.
! =============================================================================
program test_mpi_init
    use unit_test
    use mpi_communicator
    implicit none

    logical :: passed = .false.

    call mpi_comm_initialise

    if (world%size == 1) then
        print *, "MPI tests must be run with more than one process."
        call mpi_comm_finalise
        stop
    endif

    passed = (world%err == 0)

    call mpi_comm_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == 0) then
        call print_result_logical('Test MPI init', passed)
    endif

end program test_mpi_init
