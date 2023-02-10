module mpi_utils
    use mpi_communicator
    implicit none

    contains
        subroutine mpi_print(msg)
            character(*), intent(in) :: msg
            if (comm%rank == 0) then
                print *, msg
            endif
        end subroutine mpi_print

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_exit_on_error(msg)
            character(*), intent(in) :: msg
            print *, "Error on rank ", comm%rank
            print *, "    ", msg
            call MPI_Abort(comm%world, -1, comm%err)
        end subroutine mpi_exit_on_error

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_check_for_message(tag, recv_size, source)
            integer, intent(out) :: recv_size, tag, source
            type(MPI_Status)     :: status

            call MPI_probe(MPI_ANY_SOURCE,  &
                           MPI_ANY_TAG,     &
                           comm%cart,       &
                           status,          &
                           comm%err)

            source = status%MPI_SOURCE
            tag = status%MPI_TAG

            call MPI_get_count(status, MPI_DOUBLE_PRECISION, recv_size, comm%err)

        end subroutine mpi_check_for_message

end module mpi_utils
