module mpi_utils
    use mpi_communicator
    implicit none

    contains
        subroutine mpi_print(msg)
            character(*), intent(in) :: msg
            if (comm%rank == comm%master) then
                print *, msg
            endif
        end subroutine mpi_print

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_exit_on_error(msg)
            character(*), intent(in) :: msg
            print *, "Error on rank ", comm%rank
            print *, "    ", msg
            call MPI_Abort(comm%world, -1, comm%err)

            call mpi_check_for_error("in MPI_Abort of mpi_utils::mpi_exit_on_error.")
        end subroutine mpi_exit_on_error

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_check_for_message(source, tag, recv_size)
            integer, intent(in)  :: source, tag
            integer, intent(out) :: recv_size
            type(MPI_Status)     :: status

            status%MPI_SOURCE = -1
            status%MPI_TAG = -1
            status%MPI_ERROR = 0

            call MPI_probe(source,          &
                           tag,             &
                           comm%cart,       &
                           status,          &
                           comm%err)

            call mpi_check_for_error(&
                "in MPI_probe of mpi_utils::mpi_check_for_message.")

            comm%err = status%MPI_ERROR
            call mpi_check_for_error(&
                "in MPI_Status of mpi_utils::mpi_check_for_message.")

            recv_size = 0
            call MPI_get_count(status, MPI_DOUBLE_PRECISION, recv_size, comm%err)

            call mpi_check_for_error(&
                "in MPI_get_count of mpi_utils::mpi_check_for_message.")

        end subroutine mpi_check_for_message

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_check_for_error(msg)
            character(*), intent(in) :: msg
#ifndef NDEBUG
            if (.not. comm%err == MPI_SUCCESS) then
                call mpi_exit_on_error(msg)
            endif
#endif
        end subroutine mpi_check_for_error

end module mpi_utils
