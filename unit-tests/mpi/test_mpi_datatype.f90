! =============================================================================
!                       Test MPI datatypes
!
!               This unit test checks the MPI init and finalize.
! =============================================================================
program test_mpi_datatype
    use unit_test
    use mpi_environment
    use mpi_datatypes
    use mpi_collectives
    use datatypes
    implicit none

    logical :: passed = .false.
    integer(kind=int64) :: i64

    call mpi_env_initialise

    if (world%size == 1) then
        print *, "MPI tests must be run with more than one process."
        call mpi_env_finalise
        stop
    endif

    passed = (world%err == 0)


    i64 = world%rank + 1

    call mpi_blocking_reduce(i64, MPI_SUM, world)

    if (world%rank == world%root) then
        passed = (passed .and. (i64 == (world%size+1) * world%size / 2_int64))
    endif

    call mpi_env_finalise

    passed = (passed .and. (world%err == 0))

    if (world%rank == 0) then
        call print_result_logical('Test MPI datatype', passed)
    endif

end program test_mpi_datatype
