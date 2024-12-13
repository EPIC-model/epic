!==============================================================================
! This module contains user-defined MPI datatypes.
!==============================================================================
module mpi_datatypes
    ! 10 March 2024
    ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Type_create_f90_integer.3.php
    ! 13 Dec 2024
    ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node425.htm
    use mpi_f08
    use datatypes
    implicit none

    private

    type(MPI_Datatype) :: MPI_INTEGER_64BIT

    type(MPI_Datatype) :: MPI_INTEGER_LOGICAL       &
                        , MPI_INTEGER_LOGICAL_ARRAY

    public  :: mpi_datatypes_create         &
             , MPI_INTEGER_64BIT            &
             , MPI_INTEGER_LOGICAL          &
             , MPI_INTEGER_LOGICAL_ARRAY

    contains

        subroutine mpi_datatypes_create
            integer :: err

            ! 18 must be the same argument as used in selected_int_kind (in file utils/datatypes.f90)
            call MPI_Type_create_f90_integer(18, MPI_INTEGER_64BIT, err)

            call check_success(err)

            call create_integer_logical_pair

        end subroutine mpi_datatypes_create

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine  create_integer_logical_pair
            type(intlog_pair_t)            :: dtype, dtypearr(2)
            integer                        :: err
            integer                        :: blocklen(2)
            type(MPI_Datatype)             :: data_type(2)
            integer(KIND=MPI_ADDRESS_KIND) :: disp(2), base, lb, extent

            !----------------------------------------------
            ! Enable sending scalars of that pair type
            call MPI_GET_ADDRESS(dtype%ival,   disp(1), err)
            call check_success(err)

            call MPI_GET_ADDRESS(dtype%lval, disp(2), err)
            call check_success(err)

            base = disp(1)
            disp(1) = disp(1) - base
            disp(2) = disp(2) - base

            blocklen(1) = 1
            blocklen(2) = 1

            data_type(1) = MPI_INTEGER
            data_type(2) = MPI_LOGICAL

            call MPI_TYPE_CREATE_STRUCT(2, blocklen, disp, data_type, MPI_INTEGER_LOGICAL, err)
            call check_success(err)

            call MPI_TYPE_COMMIT(MPI_INTEGER_LOGICAL, err)
            call check_success(err)

            !----------------------------------------------
            ! Enable sending arrays of that pair type
            call MPI_GET_ADDRESS(dtypearr(1), disp(1), err)
            call check_success(err)

            call MPI_GET_ADDRESS(dtypearr(2), disp(2), err)
            call check_success(err)

            extent = disp(2) - disp(1)
            lb = 0
            call MPI_TYPE_CREATE_RESIZED(MPI_INTEGER_LOGICAL,       &
                                         lb,                        &
                                         extent,                    &
                                         MPI_INTEGER_LOGICAL_ARRAY, &
                                         err)

            call MPI_TYPE_COMMIT(MPI_INTEGER_LOGICAL_ARRAY, err)

        end subroutine create_integer_logical_pair

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine check_success(err)
            integer, intent(inout) :: err

            if (err /= MPI_SUCCESS) then
                print *, "Error in mpi_datatypes::mpi_datatypes_create: Unable to create user-defined MPI type."
                call MPI_Abort(MPI_COMM_WORLD, -1, err)
            endif
        end subroutine check_success

end module mpi_datatypes
