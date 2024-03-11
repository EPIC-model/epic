module mpi_datatypes
    ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_create_f90_integer.3.php
    use mpi_f08
    use datatypes
    implicit none

    type(MPI_Datatype) :: MPI_INTEGER_64BIT
    type(MPI_OP)       :: MPI_SUM_64BIT

    contains

        subroutine mpi_create_types
            integer :: err

            call MPI_Type_create_f90_integer(18, MPI_INTEGER_64BIT)

            call MPI_Op_create(my_user_function, commute=.true., op=MPI_SUM_64BIT, ierror=err)

        end subroutine mpi_create_types

        !https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node115.htm
        subroutine my_user_function( invec, inoutvec, length, dtype)
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
        end subroutine

end module mpi_datatypes
