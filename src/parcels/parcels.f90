module parcel_container
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_1d, &
                       write_h5_dataset_2d
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

        subroutine h5_write_parcels(iter)
            integer, intent(in) :: iter ! iteration

            call h5open_f(h5err)

            !
            ! write parcel data
            !

            call write_h5_dataset_2d("position", parcels%pos(1:n_parcels, :))
            call write_h5_dataset_2d("velocity", parcels%vel(1:n_parcels, :))

            if (allocated(parcels%stretch)) then
                call write_h5_dataset_1d("stretch", parcels%stretch(1:n_parcels))
            endif

            if (allocated(parcels%B11) .and. allocated(parcels%B12)) then
                call write_h5_dataset_1d("B11", parcels%B11(1:n_parcels))
                call write_h5_dataset_1d("B12", parcels%B12(1:n_parcels))
            endif

            call h5close_f(h5err)
        end subroutine h5_write_parcels

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
