! =============================================================================
! Module for common array manipulations.
! =============================================================================
module armanip
    implicit none

    private

    interface resize_array
        module procedure :: resize_array_1d
        module procedure :: resize_array_2d
    end interface resize_array

    public :: resize_array

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine resize_array_1d(attr, new_size)
            double precision, allocatable, intent(inout) :: attr(:)
            integer,           intent(in)    :: new_size
            double precision, allocatable    :: buffer(:)
            integer                          :: old_size

            old_size = size(attr)

            allocate(buffer(new_size))

            buffer(1:old_size) = attr(1:old_size)

            call move_alloc(buffer, attr)

        end subroutine resize_array_1d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine resize_array_2d(attr, new_size)
            double precision, allocatable, intent(inout) :: attr(:, :)
            integer,           intent(in)    :: new_size
            double precision, allocatable    :: buffer(:, :)
            integer                          :: shap(2)

            shap = shape(attr)

            allocate(buffer(new_size, shap(2)))

            buffer(1:shap(1), :) = attr(1:shap(1), :)

            call move_alloc(buffer, attr)

        end subroutine resize_array_2d

end module armanip
