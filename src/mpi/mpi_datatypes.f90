module mpi_datatypes
    ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_create_f90_integer.3.php
    use mpi_f08
    implicit none

    type(MPI_Datatype) :: MPI_INTEGER_64BIT

    contains

        subroutine mpi_create_types

            call MPI_Type_create_f90_integer(18, MPI_INTEGER_64BIT)

        end subroutine mpi_create_types

end module mpi_datatypes
