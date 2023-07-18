module mpi_communicator
    use mpi_f08
    use mpi_tags
    implicit none

    type :: communicator
        type(MPI_Comm) :: comm = MPI_COMM_WORLD
        integer        :: err = 0
        integer        :: rank = -1
        integer        :: size = 1
        integer        :: master = 0
    end type communicator

    type(communicator) :: world

    type :: sub_communicator
        type(MPI_Comm) :: comm = MPI_COMM_WORLD
        integer        :: err = 0
        integer        :: rank = 0
        integer        :: size = 1
    end type sub_communicator

    contains

        subroutine mpi_comm_initialise
#ifdef ENABLE_OPENMP
            integer :: provided
            call MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, world%err)

            if (.not. provided == MPI_THREAD_SERIALIZED) then
                stop
            endif
#else
            call MPI_Init(world%err)
#endif
            call MPI_Comm_size(world%comm, world%size, world%err)
            call MPI_Comm_rank(world%comm, world%rank, world%err)
            call MPI_Comm_size(world%comm, world%size, world%err)
        end subroutine mpi_comm_initialise

        subroutine mpi_comm_finalise
            logical :: flag
            call MPI_Initialized(flag, world%err)

            if (flag) then
                call MPI_Finalize(world%err)
            endif

            call MPI_Finalized(flag, world%err)

            if (.not. flag) then
                call MPI_Abort(world%comm, -1, world%err)
            endif
        end subroutine mpi_comm_finalise
end module mpi_communicator
