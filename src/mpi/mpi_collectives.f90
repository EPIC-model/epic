module mpi_collectives
    use mpi_communicator
    implicit none

    interface mpi_blocking_reduce
        module procedure :: mpi_double_reduce
        module procedure :: mpi_integer_reduce
    end interface mpi_blocking_reduce

    contains

        subroutine mpi_double_reduce(sendbuf, op)
            double precision, intent(inout) :: sendbuf(..)
            type(MPI_Op),     intent(in)    :: op

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, world%root, world%comm, world%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, world%root, world%comm, world%err)
            endif
        end subroutine mpi_double_reduce

        subroutine mpi_integer_reduce(sendbuf, op)
            integer,      intent(inout) :: sendbuf(..)
            type(MPI_Op), intent(in)    :: op

            if (world%rank == world%root) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, world%root, world%comm, world%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_INTEGER, &
                                op, world%root, world%comm, world%err)
            endif
        end subroutine mpi_integer_reduce

end module mpi_collectives
