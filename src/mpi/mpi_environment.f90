module mpi_environment
    use mpi_f08
    use mpi_tags
    use mpi_datatypes, only : mpi_datatypes_create
    use mpi_ops, only : mpi_ops_create &
                      , mpi_ops_free
    implicit none

    type :: communicator
        type(MPI_Comm) :: comm = MPI_COMM_NULL
        integer        :: err = 0
        integer        :: rank = 0
        integer        :: size = 1
        integer        :: root = 0
    end type communicator

    type(communicator) :: world

    contains

        subroutine mpi_env_initialise
#ifdef ENABLE_OPENMP
            integer :: provided
            call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, world%err)

            if (.not. provided == MPI_THREAD_FUNNELED) then
                stop
            endif
#else
            call MPI_Init(world%err)
#endif
            world%comm = MPI_COMM_WORLD

            call MPI_Comm_size(world%comm, world%size, world%err)
            call MPI_Comm_rank(world%comm, world%rank, world%err)

            call mpi_datatypes_create

            call mpi_ops_create

        end subroutine mpi_env_initialise

        subroutine mpi_env_finalise
            logical :: flag
            call MPI_Initialized(flag, world%err)

            if (flag) then
                call mpi_ops_free
                call MPI_Finalize(world%err)
            endif

            call MPI_Finalized(flag, world%err)

            if (.not. flag) then
                call MPI_Abort(world%comm, -1, world%err)
            endif
        end subroutine mpi_env_finalise
end module mpi_environment
