module mpi_communicator
    use mpi_f08
    use mpi_tags
    implicit none

    type(MPI_Comm), parameter :: comm_world = MPI_COMM_WORLD
    type(MPI_Comm)            :: comm_cart
    integer                   :: mpi_err = 0
    integer                   :: mpi_rank = 0
    integer                   :: mpi_size = 1
    integer, parameter        :: mpi_master = 0

    contains

        subroutine mpi_comm_initialise
            call MPI_Init(mpi_err)
            call MPI_Comm_size(comm_world, mpi_size, mpi_err)
            call MPI_Comm_rank(comm_world, mpi_rank, mpi_err)
            call MPI_Comm_size(comm_world, mpi_size, mpi_err)
        end subroutine mpi_comm_initialise

        subroutine mpi_comm_finalise
            logical :: flag
            call MPI_Initialized(flag, mpi_err)

            if (flag) then
                call MPI_Finalize(mpi_err)
            endif

            call MPI_Finalized(flag, mpi_err)

            if (.not. flag) then
                call MPI_Abort(comm_world, -1, mpi_err)
            endif
        end subroutine mpi_comm_finalise
end module mpi_communicator
