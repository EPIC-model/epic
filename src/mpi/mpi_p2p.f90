module mpi_p2p
    use mpi_communicator
    implicit none

    interface mpi_non_blocking_send
        module procedure :: mpi_double_isend
        module procedure :: mpi_integer_isend
    end interface mpi_non_blocking_send

    interface mpi_non_blocking_recv
        module procedure :: mpi_double_irecv
        module procedure :: mpi_integer_irecv
    end interface mpi_non_blocking_recv

    interface mpi_blocking_send
        module procedure :: mpi_double_send
        module procedure :: mpi_integer_send
    end interface mpi_blocking_send

    interface mpi_blocking_recv
        module procedure :: mpi_double_recv
        module procedure :: mpi_integer_recv
    end interface mpi_blocking_recv


    contains
        !
        ! non-blocking send/receive
        !
        ! Uses 'asynchronous', see e.g. https://www.open-mpi.org/doc/current/man3/MPI_Isend.3.php
        ! and https://stackoverflow.com/questions/19455051/using-mpi-send-recv-to-handle-chunk-of-multi-dim-array-in-fortran-90 (23 November 2021)
        !
        subroutine mpi_double_isend(data, dest, tag)
            double precision, intent(in), asynchronous :: data(..)
            integer,          intent(in)               :: dest
            integer,          intent(in)               :: tag
            type(MPI_Request)                          :: request

            call MPI_Isend(data(1:size(data)),      &
                           size(data),              &
                           MPI_DOUBLE_PRECISION,    &
                           dest,                    &
                           tag,                     &
                           comm%world,              &
                           request,                 &
                           comm%err)
        end subroutine mpi_double_isend

        subroutine mpi_integer_isend(data, dest, tag)
            integer, intent(in), asynchronous :: data(..)
            integer, intent(in)               :: dest
            integer, intent(in)               :: tag
            type(MPI_Request)                 :: request

            call MPI_Isend(data(1:size(data)),  &
                           size(data),          &
                           MPI_INTEGER,         &
                           dest,                &
                           tag,                 &
                           comm%world,          &
                           request,             &
                           comm%err)
        end subroutine mpi_integer_isend


        subroutine mpi_double_irecv(data, source, tag)
            double precision, intent(out), asynchronous :: data(..)
            integer,          intent(in)                :: source
            integer,          intent(in)                :: tag
            type(MPI_Request)                           :: request

            call MPI_Irecv(data(1:size(data)),      &
                           size(data),              &
                           MPI_DOUBLE_PRECISION,    &
                           source,                  &
                           tag,                     &
                           comm%world,              &
                           request,                 &
                           comm%err)
        end subroutine mpi_double_irecv

        subroutine mpi_integer_irecv(data, source, tag)
            integer, intent(out), asynchronous :: data(..)
            integer, intent(in)                :: source
            integer, intent(in)                :: tag
            type(MPI_Request)                  :: request

            call MPI_Irecv(data(1:size(data)),  &
                           size(data),          &
                           MPI_INTEGER,         &
                           source,              &
                           tag,                 &
                           comm%world,          &
                           request,             &
                           comm%err)
        end subroutine mpi_integer_irecv


        !
        ! blocking send/receive
        !
        subroutine mpi_double_send(data, dest, tag)
            double precision, intent(in) :: data(..)
            integer,          intent(in) :: dest
            integer,          intent(in) :: tag

            call MPI_Send(data(1:size(data)),   &
                          size(data),           &
                          MPI_DOUBLE_PRECISION, &
                          dest,                 &
                          tag,                  &
                          comm%world,           &
                          comm%err)
        end subroutine mpi_double_send

        subroutine mpi_integer_send(data, dest, tag)
            integer, intent(in) :: data(..)
            integer, intent(in) :: dest
            integer, intent(in) :: tag

            call MPI_Send(data(1:size(data)),   &
                          size(data),           &
                          MPI_INTEGER,          &
                          dest,                 &
                          tag,                  &
                          comm%world,           &
                          comm%err)
        end subroutine mpi_integer_send


        subroutine mpi_double_recv(data, source, tag)
            double precision, intent(out) :: data(..)
            integer,          intent(in)  :: source
            integer,          intent(in)  :: tag
            type(MPI_Status)              :: status

            call MPI_Recv(data(1:size(data)),   &
                          size(data),           &
                          MPI_DOUBLE_PRECISION, &
                          source,               &
                          tag,                  &
                          comm%world,           &
                          status,               &
                          comm%err)
        end subroutine mpi_double_recv

        subroutine mpi_integer_recv(data, source, tag)
            integer, intent(out) :: data(..)
            integer, intent(in)  :: source
            integer, intent(in)  :: tag
            type(MPI_Status)     :: status

            call MPI_Recv(data(1:size(data)),   &
                          size(data),           &
                          MPI_INTEGER,          &
                          source,               &
                          tag,                  &
                          comm%world,           &
                          status,               &
                          comm%err)
        end subroutine mpi_integer_recv

end module mpi_p2p
