module mpi_communicator
    use mpi_f08
    use mpi_tags
    implicit none

    type :: parallel_communicator
        type(MPI_Comm) :: world = MPI_COMM_WORLD
        type(MPI_Comm) :: cart
        integer        :: err = 0
        integer        :: rank = 0
        integer        :: size = 1
        integer        :: master = 0
    end type parallel_communicator

    type(parallel_communicator) :: comm

    type :: sub_communicator(maxdims)
        integer, len                :: maxdims ! length of vector coord
        integer, dimension(maxdims) :: coord
        type(MPI_Comm)              :: comm = MPI_COMM_WORLD
        integer                     :: err = 0
        integer                     :: rank = 0
        integer                     :: size = 1
    end type sub_communicator

    contains

        subroutine mpi_comm_initialise
            call MPI_Init(comm%err)
            call MPI_Comm_size(comm%world, comm%size, comm%err)
            call MPI_Comm_rank(comm%world, comm%rank, comm%err)
            call MPI_Comm_size(comm%world, comm%size, comm%err)
        end subroutine mpi_comm_initialise

        subroutine mpi_comm_finalise
            logical :: flag
            call MPI_Initialized(flag, comm%err)

            if (flag) then
                call MPI_Finalize(comm%err)
            endif

            call MPI_Finalized(flag, comm%err)

            if (.not. flag) then
                call MPI_Abort(comm%world, -1, comm%err)
            endif
        end subroutine mpi_comm_finalise
end module mpi_communicator
