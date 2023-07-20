module mpi_collectives
    use mpi_environment
    implicit none

    interface mpi_blocking_reduce
        module procedure :: mpi_double_reduce
        module procedure :: mpi_integer_reduce
    end interface mpi_blocking_reduce

    contains

        subroutine mpi_double_reduce(sendbuf, op, comm)
            double precision,   intent(inout) :: sendbuf(..)
            type(MPI_Op),       intent(in)    :: op
            type(communicator), intent(inout) :: comm

            if (comm%rank == comm%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, comm%root, comm%comm, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, comm%root, comm%comm, comm%err)
            endif
        end subroutine mpi_double_reduce

        subroutine mpi_integer_reduce(sendbuf, op, comm)
            integer,      intent(inout) :: sendbuf(..)
            type(MPI_Op), intent(in)    :: op

            if (comm%rank == comm%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, comm%root, comm%comm, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, comm%root, comm%comm, comm%err)
            endif
        end subroutine mpi_integer_reduce

end module mpi_collectives
