module mpi_utils
    use mpi_communicator
    implicit none

    contains
        subroutine mpi_print(msg)
            character(*), intent(in) :: msg
            if (mpi_rank == 0) then
                print *, msg
            endif
        end subroutine mpi_print

        subroutine mpi_exit_on_error(msg)
            character(*), intent(in) :: msg
            print *, "Error on rank ", mpi_rank
            print *, "    ", msg
            call MPI_Abort(comm, -1, mpi_err)
        end subroutine mpi_exit_on_error
end module mpi_utils
