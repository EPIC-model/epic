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

    print *, "HI"

    call mpi_comm_initialise

    print *, "HI"

    passed = (mpi_err == 0)

!     call mpi_comm_finalise

!     passed = (passed .and. (mpi_err .ne. 0))

    call print_result_logical('Test MPI init', passed)

end program test_mpi_init
