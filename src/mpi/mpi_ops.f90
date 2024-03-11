!==============================================================================
! This module contains user-defined MPI operators for reductions
! (MPI_Reduce, MPI_Allreduce, MPI_Reduce_scatter and MPI_Scan).
! See also https://www.open-mpi.org/doc/v3.1/man3/MPI_Op_create.3.php
!==============================================================================
module mpi_ops
    use mpi_f08
    use datatypes, only : int64
#ifndef NDEBUG
    use mpi_datatypes, only : MPI_INTEGER_64BIT
#endif
    implicit none

    type(MPI_Op) :: MPI_SUM_64BIT = MPI_OP_NULL

    private :: mpi_op_sum_integer_64bit

    contains

        subroutine mpi_ops_create
            integer :: err

            call MPI_Op_create(user_fn=mpi_op_sum_integer_64bit,    &
                               commute=.true.,                      &
                               op=MPI_SUM_64BIT,                    &
                               ierror=err)

            if (err /= MPI_SUCCESS) then
                print *, "Error in mpi_ops::mpi_ops_create: Unable to create user-defined MPI operator."
                call MPI_Abort(MPI_COMM_WORLD, -1, err)
            endif

        end subroutine mpi_ops_create

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine mpi_ops_free
            integer :: err

            call MPI_Op_free(MPI_SUM_64BIT, err)

            if (err /= MPI_SUCCESS) then
                print *, "Error in mpi_ops::mpi_ops_free: Unable to free user-defined MPI operator."
                call MPI_Abort(MPI_COMM_WORLD, -1, err)
            endif

        end subroutine mpi_ops_free

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! 11 March 2024
        !https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node115.htm
        subroutine mpi_op_sum_integer_64bit(invec, inoutvec, length, dtype)
            use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
            type(c_ptr), value           :: invec, inoutvec
            integer                      :: length
            type(MPI_Datatype)           :: dtype
            integer(kind=int64), pointer :: invec64bit(:), inoutvec64bit(:)
#ifndef NDEBUG
            if (dtype /= MPI_INTEGER_64BIT) then
                call MPI_Abort(MPI_COMM_WORLD, -1, length)
            endif
#endif
                call c_f_pointer(invec, invec64bit, (/ length /) )
                call c_f_pointer(inoutvec, inoutvec64bit, (/ length /) )
                inoutvec64bit = invec64bit + inoutvec64bit
        end subroutine mpi_op_sum_integer_64bit

end module mpi_ops
