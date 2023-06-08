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

        subroutine resize_array_1d(attr, new_size, copy_size)
            double precision, allocatable, intent(inout) :: attr(:)
            integer,                       intent(in)    :: new_size
            integer, optional,             intent(in)    :: copy_size
            double precision, allocatable                :: buffer(:)
            integer                                      :: old_size

            old_size = size(attr)
            if (present(copy_size)) then
                old_size = copy_size
            endif

            if (new_size < old_size) then
                call mpi_exit_on_error("armanip::resize_array_1d: new_size < old_size")
            endif

            allocate(buffer(new_size))

            buffer(1:old_size) = attr(1:old_size)

            call move_alloc(buffer, attr)

        end subroutine resize_array_1d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine resize_array_2d(attr, new_size, copy_size)
            double precision, allocatable, intent(inout) :: attr(:, :)
            integer,                       intent(in)    :: new_size
            integer, optional,             intent(in)    :: copy_size
            double precision, allocatable                :: buffer(:, :)
            integer                                      :: ncomp, old_size

            ncomp = size(attr, dim=1)
            old_size = size(attr, dim=2)

            if (present(copy_size)) then
                old_size = copy_size
            endif

            if (new_size < old_size) then
                call mpi_exit_on_error("armanip::resize_array_2d: new_size < old_size")
            endif

            allocate(buffer(ncomp, new_size))

            buffer(:, 1:old_size) = attr(:, 1:old_size)

            call move_alloc(buffer, attr)

        end subroutine resize_array_2d

end module armanip
