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
!             call MPI_Ireduce(sendbuf, recvbuf, size(recvbuf), MPI_DOUBLE_PRECISION, &
!                              op, 0, comm%world, request, comm%err)
!         end subroutine mpi_double_ireduce

!         subroutine mpi_integer_ireduce(sendbuf, recvbuf, op)
!             integer,           intent(in),  asynchronous :: sendbuf(..)
!             integer,           intent(out), asynchronous :: recvbuf(..)
!             type(MPI_Op),      intent(in)                :: op
!             type(MPI_Request)                            :: request
!             call MPI_Ireduce(sendbuf, recvbuf, size(recvbuf), MPI_INT, &
!                              op, 0, comm%world, request, comm%err)
!         end subroutine mpi_integer_ireduce

        subroutine mpi_double_reduce(sendbuf, op)
            double precision, intent(inout) :: sendbuf(..)
            type(MPI_Op),     intent(in)    :: op

            if (comm%rank == comm%master) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, comm%master, comm%world, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_DOUBLE_PRECISION, &
                                op, comm%master, comm%world, comm%err)
            endif
        end subroutine mpi_double_reduce

        subroutine mpi_integer_reduce(sendbuf, op)
            integer,      intent(inout) :: sendbuf(..)
            type(MPI_Op), intent(in)    :: op

            if (comm%rank == comm%master) then
                call MPI_Reduce(MPI_IN_PLACE, sendbuf, size(sendbuf), MPI_INT, &
                                op, comm%master, comm%world, comm%err)
            else
                call MPI_Reduce(sendbuf, sendbuf, size(sendbuf), MPI_INT, &
                                op, comm%master, comm%world, comm%err)
            endif
        end subroutine mpi_integer_reduce

end module mpi_collectives
