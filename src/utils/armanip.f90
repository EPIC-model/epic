! =============================================================================
! Module for common array manipulations.
! =============================================================================
module armanip
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    private

    interface resize_array
        module procedure :: resize_array_1d
        module procedure :: resize_array_2d
    end interface resize_array

    public :: resize_array

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Resize a one-dimensional array
        ! @param[inout] attr the array
        ! @param[in] new_size of the array
        ! @param[in] n_copy number of elements to copy over
        subroutine resize_array_1d(attr, new_size, n_copy)
            double precision, allocatable, intent(inout) :: attr(:)
            integer,                       intent(in)    :: new_size
            integer, optional,             intent(in)    :: n_copy
            double precision, allocatable                :: buffer(:)
            integer                                      :: copy_size

            copy_size = size(attr)
            if (present(n_copy)) then
                copy_size = n_copy
            endif

            if (new_size < copy_size) then
                call mpi_exit_on_error("armanip::resize_array_1d: new_size < copy_size")
            endif

            allocate(buffer(new_size))

            buffer(1:copy_size) = attr(1:copy_size)

            call move_alloc(buffer, attr)

        end subroutine resize_array_1d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Resize a two-dimensional array
        ! @param[inout] attr the array
        ! @param[in] new_size of the array
        ! @param[in] n_copy number of elements to copy over
        subroutine resize_array_2d(attr, new_size, n_copy)
            double precision, allocatable, intent(inout) :: attr(:, :)
            integer,                       intent(in)    :: new_size
            integer, optional,             intent(in)    :: n_copy
            double precision, allocatable                :: buffer(:, :)
            integer                                      :: ncomp, copy_size

            ncomp = size(attr, dim=1)
            copy_size = size(attr, dim=2)

            if (present(n_copy)) then
                copy_size = n_copy
            endif

            if (new_size < copy_size) then
                call mpi_exit_on_error("armanip::resize_array_2d: new_size < copy_size")
            endif

            allocate(buffer(ncomp, new_size))

            buffer(:, 1:copy_size) = attr(:, 1:copy_size)

            call move_alloc(buffer, attr)

        end subroutine resize_array_2d

end module armanip
