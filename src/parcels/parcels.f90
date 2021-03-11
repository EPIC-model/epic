module parcel_container
    implicit none

    integer :: n_parcels

    type attribute_container_type
        double precision, allocatable, dimension(:) :: &
            stretch,    &
            B11, B12       ! B matrix entries

        double precision, allocatable, dimension(:, :) :: &
            pos,        & ! positions
            vel           ! velocitues
    end type attribute_container_type

    type(attribute_container_type) parcels


    contains

        subroutine split(threshold)
            double precision, intent(in) :: threshold


        end subroutine split

        subroutine alloc_parcel_mem(num)
            integer, intent(in) :: num

            allocate(parcels%pos(num, 2))
            allocate(parcels%vel(num, 2))
            allocate(parcels%stretch(num))
            allocate(parcels%B11(num))
            allocate(parcels%B12(num))
        end subroutine alloc_parcel_mem

        subroutine dealloc_parcel_mem
            deallocate(parcels%pos)
            deallocate(parcels%vel)
            deallocate(parcels%stretch)
            deallocate(parcels%B11)
            deallocate(parcels%B12)
        end subroutine dealloc_parcel_mem

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcel_container
