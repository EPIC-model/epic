module parcels
    implicit none

    integer :: n_parcels

    double precision, allocatable, dimension(:) :: &
        stretch,    &
        B11, B12       ! B matrix entries

    double precision, allocatable, dimension(:, :) :: &
        pos,        & ! positions
        vel           ! velocitues


    contains

        subroutine split(threshold)
            double precision, intent(in) :: threshold


        end subroutine split

        subroutine alloc_parcel_mem(num)
            integer, intent(in) :: num

            allocate(pos(num, 2))
            allocate(vel(num, 2))
            allocate(stretch(num))
            allocate(B11(num))
            allocate(B12(num))
        end subroutine alloc_parcel_mem

        subroutine dealloc_parcel_mem
            deallocate(pos)
            deallocate(vel)
            deallocate(stretch)
            deallocate(B11)
            deallocate(B12)
        end subroutine dealloc_parcel_mem

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcels
