module mpi_collectives
    use mpi_communicator
    implicit none

!     interface mpi_non_blocking_reduce
!         module procedure :: mpi_double_ireduce
!         module procedure :: mpi_integer_ireduce
!     end interface mpi_non_blocking_reduce

    interface mpi_blocking_reduce
        module procedure :: mpi_double_reduce
        module procedure :: mpi_integer_reduce
    end interface mpi_blocking_reduce

    contains
!         subroutine mpi_double_ireduce(sendbuf, recvbuf, op)
!             double precision,  intent(in),  asynchronous :: sendbuf(..)
!             double precision,  intent(out), asynchronous :: recvbuf(..)
!             type(MPI_Op),      intent(in)                :: op
!             type(MPI_Request)                            :: request
!             call MPI_Ireduce(sendbuf, recvbuf, size(recvbuf), MPI_DOUBLE, &
!                              op, 0, comm_world, request, mpi_err)
!         end subroutine mpi_double_ireduce

!         subroutine mpi_integer_ireduce(sendbuf, recvbuf, op)
!             integer,           intent(in),  asynchronous :: sendbuf(..)
!             integer,           intent(out), asynchronous :: recvbuf(..)
!             type(MPI_Op),      intent(in)                :: op
!             type(MPI_Request)                            :: request
!             call MPI_Ireduce(sendbuf, recvbuf, size(recvbuf), MPI_INT, &
!                              op, 0, comm_world, request, mpi_err)
!         end subroutine mpi_integer_ireduce

        subroutine mpi_double_reduce(sendbuf, recvbuf, op)
            double precision, intent(in)  :: sendbuf(..)
            double precision, intent(out) :: recvbuf(..)
            type(MPI_Op),     intent(in)  :: op
            call MPI_Reduce(sendbuf, recvbuf, size(recvbuf), MPI_DOUBLE, &
                            op, 0, comm_world, mpi_err)
        end subroutine mpi_double_reduce

        subroutine mpi_integer_reduce(sendbuf, recvbuf, op)
            integer,      intent(in)  :: sendbuf(..)
            integer,      intent(out) :: recvbuf(..)
            type(MPI_Op), intent(in)  :: op
            call MPI_Reduce(sendbuf, recvbuf, size(recvbuf), MPI_INT, &
                            op, 0, comm_world, mpi_err)
        end subroutine mpi_integer_reduce

end module mpi_collectives
