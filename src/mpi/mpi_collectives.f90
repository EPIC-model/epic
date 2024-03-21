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

    interface mpi_neighbor_allreduce
        module procedure :: mpi_neighbor_allreduce_logical
        module procedure :: mpi_neighbor_allreduce_integer
    end interface mpi_neighbor_allreduce

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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_neighbor_allreduce_logical(sendbuf, recvbuf, op)
            logical,            intent(in)    :: sendbuf
            logical,            intent(inout) :: recvbuf
            type(MPI_Op),       intent(in)    :: op
            type(MPI_Request)                 :: requests(16)
            logical                           :: buf(8)
            integer                           :: n

            do n = 1, 8
                call MPI_Isend(sendbuf,                 &
                               1,                       &
                               MPI_LOGICAL,             &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of mpi_neighbor_allreduce_logical.")

                call MPI_Irecv(buf(n),                  &
                              1,                        &
                              MPI_LOGICAL,              &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              requests(8+n),            &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Irecv of mpi_neighbor_allreduce_logical.")

            enddo

            call MPI_Waitall(16,                  &
                             requests,            &
                             MPI_STATUSES_IGNORE, &
                             cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of mpi_neighbor_allreduce_logical.")

            ! Combine locally:
            recvbuf = sendbuf
            do n = 1, 8
                call MPI_Reduce_local(buf(n), recvbuf, 1, MPI_LOGICAL, op, cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Reduce_local of mpi_neighbor_allreduce_logical.")
            enddo

        end subroutine mpi_neighbor_allreduce_logical

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_neighbor_allreduce_integer(sendbuf, recvbuf, op)
            integer,            intent(in)    :: sendbuf
            integer,            intent(inout) :: recvbuf
            type(MPI_Op),       intent(in)    :: op
            type(MPI_Request)                 :: requests(16)
            integer                           :: buf(8)
            integer                           :: n

            do n = 1, 8
                call MPI_Isend(sendbuf,                 &
                               1,                       &
                               MPI_INTEGER,             &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of mpi_neighbor_allreduce_integer.")

                call MPI_Irecv(buf(n),                  &
                              1,                        &
                              MPI_INTEGER,              &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              requests(8+n),            &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Irecv of mpi_neighbor_allreduce_integer.")

            enddo

            call MPI_Waitall(16,                  &
                             requests,            &
                             MPI_STATUSES_IGNORE, &
                             cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of mpi_neighbor_allreduce_integer.")

            ! Combine locally:
            recvbuf = sendbuf
            do n = 1, 8
                call MPI_Reduce_local(buf(n), recvbuf, 1, MPI_INTEGER, op, cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Reduce_local of mpi_neighbor_allreduce_integer.")
            enddo

        end subroutine mpi_neighbor_allreduce_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_neighbor_barrier
            logical sendbuf, recvbuf
            sendbuf = .true.
            recvbuf = .false.
            call mpi_neighbor_allreduce_logical(sendbuf, recvbuf, MPI_LAND)

            if (.not. recvbuf) then
                call mpi_check_for_error(cart, &
                    "in mpi_neighbor_barrier: Barrier did not work.")
            endif

        end subroutine mpi_neighbor_barrier


end module mpi_collectives
