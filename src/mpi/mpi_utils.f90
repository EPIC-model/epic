module mpi_utils
    use mpi_environment
    implicit none

    contains
        subroutine mpi_print(msg)
            character(*), intent(in) :: msg
            if (world%rank == world%root) then
                print *, msg
            endif
        end subroutine mpi_print

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_exit_on_error(msg, rank)
            character(*), optional, intent(in) :: msg
            integer,      optional, intent(in) :: rank
            integer                            :: r

            if (present(msg)) then
                r = world%rank
                if (present(rank)) then
                    r = rank
                endif
                print *, "Error on rank ", r, msg
            endif
            call MPI_Abort(world%comm, -1, world%err)

            call mpi_check_for_error(world, &
                "in MPI_Abort of mpi_utils::mpi_exit_on_error.")
        end subroutine mpi_exit_on_error

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_stop(msg)
            character(*), optional, intent(in) :: msg
            if (present(msg)) then
                call mpi_print(msg)
            endif
            call mpi_env_finalise
            stop
        end subroutine mpi_stop

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_check_for_message(source, tag, recv_size, comm)
            integer,            intent(in)    :: source, tag
            integer,            intent(out)   :: recv_size
            type(communicator), intent(inout) :: comm
            type(MPI_Status)                  :: status

            status%MPI_SOURCE = -1
            status%MPI_TAG = -1
            status%MPI_ERROR = 0

            call MPI_probe(source,          &
                           tag,             &
                           comm%comm,       &
                           status,          &
                           comm%err)

            call mpi_check_for_error(comm, &
                "in MPI_probe of mpi_utils::mpi_check_for_message.")

            comm%err = status%MPI_ERROR
            call mpi_check_for_error(comm, &
                "in MPI_Status of mpi_utils::mpi_check_for_message.")

            recv_size = 0
            call MPI_get_count(status, MPI_DOUBLE_PRECISION, recv_size, comm%err)

            call mpi_check_for_error(comm, &
                "in MPI_get_count of mpi_utils::mpi_check_for_message.")

        end subroutine mpi_check_for_message

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_check_for_error(comm, msg)
            type(communicator), intent(in) :: comm
            character(*),       intent(in) :: msg
#ifndef NDEBUG
            if (.not. comm%err == MPI_SUCCESS) then
                call mpi_exit_on_error(msg, comm%rank)
            endif
#endif
        end subroutine mpi_check_for_error

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Check the memory model for the MPI RMA communication and stop if
        ! not MPI_WIN_UNIFIED. There exist two memory models:
        !  - MPI_WIN_SEPARATE --> public copy /= private copy
        !  - MPI_WIN_UNIFIED  --> public copy == private copy
        !
        ! EPIC only supports MPI_WIN_UNIFIED at the moment.
        ! In case of MPI_WIN_SEPARATE, the public and private copy must be
        ! synchronised with MPI_Win_sync. Every time the private copy changes,
        ! it must be made visible with MPI_WIN_SYNC. A call to MPI_Barrier may/must be
        ! done afterwards. The call to MPI_Win_sync must be within an epoch, e.g.
        !       call MPI_lock_all(win, ...)
        !       call MPI_Win_sync(win, ...)
        !       call MPI_Barrier(...)
        !       call MPI_unlock_all(win, ...)
        subroutine mpi_check_rma_window_model(win)
            type(MPI_Win), intent(in)      :: win
            integer(kind=MPI_ADDRESS_KIND) :: memory_model
            logical                        :: flag

            ! pre-set
            flag = .false.
            memory_model = MPI_WIN_SEPARATE

            ! check type of window: MPI_WIN_SEPARATE or MPI_WIN_UNIFIED
            call MPI_Win_get_attr(win,              &
                                  MPI_WIN_MODEL,    &
                                  memory_model,     &
                                  flag,             &
                                  world%err);

            if (.not. flag) then
                call mpi_exit_on_error(&
                    "in mpi_check_rma_window_model: No attribute is associated with the key.")
            endif

            if (memory_model /= MPI_WIN_UNIFIED) then
                call mpi_stop("No unified window memory model.")
            endif
        end subroutine mpi_check_rma_window_model

end module mpi_utils
