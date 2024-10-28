!==============================================================================
! This module contains user-defined MPI datatypes.
!==============================================================================
module mpi_datatypes
    ! 10 March 2024
    ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_create_f90_integer.3.php
    use mpi_f08
    implicit none

    type(MPI_Datatype) :: MPI_INTEGER_64BIT

    contains

        subroutine mpi_datatypes_create
            integer :: err

            ! 18 must be the same argument as used in selected_int_kind (in file utils/datatypes.f90)
            call MPI_Type_create_f90_integer(18, MPI_INTEGER_64BIT, err)

            if (err /= MPI_SUCCESS) then
                print *, "Error in mpi_datatypes::mpi_datatypes_create: Unable to create user-defined MPI type."
                call MPI_Abort(MPI_COMM_WORLD, -1, err)
            endif

        end subroutine mpi_datatypes_create

end module mpi_datatypes
