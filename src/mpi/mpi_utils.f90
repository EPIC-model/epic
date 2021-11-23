module mpi_utils
    use mpi_communicator
    implicit none

    interface mpi_non_blocking_send
        module procedure :: mpi_dscalar_isend
        module procedure :: mpi_iscalar_isend
        module procedure :: mpi_darray_isend
        module procedure :: mpi_iarray_isend
    end interface mpi_non_blocking_send

    interface mpi_non_blocking_recv
        module procedure :: mpi_dscalar_irecv
        module procedure :: mpi_iscalar_irecv
        module procedure :: mpi_darray_irecv
        module procedure :: mpi_iarray_irecv
    end interface mpi_non_blocking_recv

    interface mpi_blocking_send
        module procedure :: mpi_dscalar_send
        module procedure :: mpi_iscalar_send
        module procedure :: mpi_darray_send
        module procedure :: mpi_iarray_send
    end interface mpi_blocking_send

    interface mpi_blocking_recv
        module procedure :: mpi_dscalar_recv
        module procedure :: mpi_iscalar_recv
        module procedure :: mpi_darray_recv
        module procedure :: mpi_iarray_recv
    end interface mpi_blocking_recv


    contains
        !
        ! non-blocking send/receive
        !
        ! Uses 'asynchronous', see e.g. https://www.open-mpi.org/doc/current/man3/MPI_Isend.3.php
        ! and https://stackoverflow.com/questions/19455051/using-mpi-send-recv-to-handle-chunk-of-multi-dim-array-in-fortran-90 (23 November 2021)
        !
        subroutine mpi_dscalar_isend(data, dest, tag)
            double precision, intent(in), asynchronous :: data
            integer,          intent(in)               :: dest
            integer,          intent(in)               :: tag
            type(MPI_Request)                          :: request
            call MPI_Isend(data, 1, MPI_DOUBLE, dest, tag, comm, request, mpi_err)
        end subroutine mpi_dscalar_isend

        subroutine mpi_iscalar_isend(data, dest, tag)
            integer, intent(in), asynchronous :: data
            integer, intent(in)               :: dest
            integer, intent(in)               :: tag
            type(MPI_Request)                 :: request
            call MPI_Isend(data, 1, MPI_INT, dest, tag, comm, request, mpi_err)
        end subroutine mpi_iscalar_isend

        subroutine mpi_darray_isend(data, dest, tag)
            double precision, intent(in), asynchronous :: data(:)
            integer,          intent(in)               :: dest
            integer,          intent(in)               :: tag
            type(MPI_Request)                          :: request
            call MPI_Isend(data, size(data), MPI_DOUBLE, dest, tag, comm, request, mpi_err)
        end subroutine mpi_darray_isend

        subroutine mpi_iarray_isend(data, dest, tag)
            integer, intent(in), asynchronous :: data(:)
            integer, intent(in)               :: dest
            integer, intent(in)               :: tag
            type(MPI_Request)                 :: request
            call MPI_Isend(data, size(data), MPI_INT, dest, tag, comm, request, mpi_err)
        end subroutine mpi_iarray_isend



        subroutine mpi_dscalar_irecv(data, source, tag)
            double precision, intent(out), asynchronous :: data
            integer,          intent(in)                :: source
            integer,          intent(in)                :: tag
            type(MPI_Request)                           :: request
            call MPI_Irecv(data, 1, MPI_DOUBLE, source, tag, comm, request, mpi_err)
        end subroutine mpi_dscalar_irecv

        subroutine mpi_iscalar_irecv(data, source, tag)
            integer, intent(out), asynchronous :: data
            integer, intent(in)                :: source
            integer, intent(in)                :: tag
            type(MPI_Request)                  :: request
            call MPI_Irecv(data, 1, MPI_INT, source, tag, comm, request, mpi_err)
        end subroutine mpi_iscalar_irecv

        subroutine mpi_darray_irecv(data, source, tag)
            double precision, intent(out), asynchronous :: data(:)
            integer,          intent(in)                :: source
            integer,          intent(in)                :: tag
            type(MPI_Request)                           :: request
            call MPI_Irecv(data, size(data), MPI_DOUBLE, source, tag, comm, request, mpi_err)
        end subroutine mpi_darray_irecv

        subroutine mpi_iarray_irecv(data, source, tag)
            integer, intent(out), asynchronous :: data(:)
            integer, intent(in)                :: source
            integer, intent(in)                :: tag
            type(MPI_Request)                  :: request
            call MPI_Irecv(data, size(data), MPI_INT, source, tag, comm, request, mpi_err)
        end subroutine mpi_iarray_irecv


        !
        ! blocking send/receive
        !
        subroutine mpi_dscalar_send(data, dest, tag)
            double precision, intent(in) :: data
            integer,          intent(in) :: dest
            integer,          intent(in) :: tag
            call MPI_Send(data, 1, MPI_DOUBLE, dest, tag, comm, mpi_err)
        end subroutine mpi_dscalar_send

        subroutine mpi_iscalar_send(data, dest, tag)
            integer, intent(in) :: data
            integer, intent(in) :: dest
            integer, intent(in) :: tag
            call MPI_Send(data, 1, MPI_INT, dest, tag, comm, mpi_err)
        end subroutine mpi_iscalar_send

        subroutine mpi_darray_send(data, dest, tag)
            double precision, intent(in) :: data(:)
            integer,          intent(in) :: dest
            integer,          intent(in) :: tag
            call MPI_Send(data, size(data), MPI_DOUBLE, dest, tag, comm, mpi_err)
        end subroutine mpi_darray_send

        subroutine mpi_iarray_send(data, dest, tag)
            integer, intent(in) :: data(:)
            integer, intent(in) :: dest
            integer, intent(in) :: tag
            call MPI_Send(data, size(data), MPI_INT, dest, tag, comm, mpi_err)
        end subroutine mpi_iarray_send



        subroutine mpi_dscalar_recv(data, source, tag)
            double precision, intent(out) :: data
            integer,          intent(in)  :: source
            integer,          intent(in)  :: tag
            type(MPI_Status)              :: status
            call MPI_Recv(data, 1, MPI_DOUBLE, source, tag, comm, status, mpi_err)
        end subroutine mpi_dscalar_recv

        subroutine mpi_iscalar_recv(data, source, tag)
            integer, intent(out) :: data
            integer, intent(in)  :: source
            integer, intent(in)  :: tag
            type(MPI_Status)     :: status
            call MPI_Recv(data, 1, MPI_INT, source, tag, comm, status, mpi_err)
        end subroutine mpi_iscalar_recv

        subroutine mpi_darray_recv(data, source, tag)
            double precision, intent(out) :: data(:)
            integer,          intent(in)  :: source
            integer,          intent(in)  :: tag
            type(MPI_Status)              :: status
            call MPI_Recv(data, size(data), MPI_DOUBLE, source, tag, comm, status, mpi_err)
        end subroutine mpi_darray_recv

        subroutine mpi_iarray_recv(data, source, tag)
            integer, intent(out) :: data(:)
            integer, intent(in)  :: source
            integer, intent(in)  :: tag
            type(MPI_Status)     :: status
            call MPI_Recv(data, size(data), MPI_INT, source, tag, comm, status, mpi_err)
        end subroutine mpi_iarray_recv

end module mpi_utils
