module mpi_collectives
    use datatypes, only : int64
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_environment
    use mpi_utils, only : mpi_stop, mpi_check_for_error
    use mpi_layout
    implicit none

    interface mpi_blocking_reduce
        module procedure :: mpi_double_reduce
        module procedure :: mpi_integer_reduce
        module procedure :: mpi_integer64bit_reduce
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_integer_reduce(sendbuf, op, comm)
            integer,            intent(inout) :: sendbuf(..)
            type(MPI_Op),       intent(in)    :: op
            type(communicator), intent(inout) :: comm

            if (comm%rank == comm%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, comm%root, comm%comm, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, comm%root, comm%comm, comm%err)
            endif
        end subroutine mpi_integer_reduce

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_integer64bit_reduce(sendbuf, op, comm)
            integer(kind=int64), intent(inout) :: sendbuf(..)
            type(MPI_Op),        intent(in)    :: op
            type(communicator),  intent(inout) :: comm

            if (op /= MPI_SUM) then
                call mpi_stop("Only MPI_SUM supported for 64-bit integers!")
            endif

            if (comm%rank == comm%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_INTEGER_64BIT, &
                                MPI_SUM_64BIT, comm%root, comm%comm, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_INTEGER_64BIT, &
                                MPI_SUM_64BIT, comm%root, comm%comm, comm%err)
            endif
        end subroutine mpi_integer64bit_reduce

end module mpi_collectives
