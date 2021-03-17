module parcel_container
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_1d, &
                       write_h5_dataset_2d, &
                       open_h5_group,       &
                       get_step_group_name
    implicit none

    integer :: n_parcels

    type parcel_container_type
        double precision, allocatable, dimension(:) :: &
            stretch,    &
            volume

        double precision, allocatable, dimension(:, :) :: &
            position,   &
            velocity,   &
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        subroutine write_h5_parcels(iter)
            integer, intent(in)        :: iter ! iteration
            integer(hid_t)             :: group
            integer(hid_t)             :: step_group
            character(:), allocatable  :: step
            character(:), allocatable  :: name
            logical                    :: link_exists = .false.

            step = trim(get_step_group_name(iter))

            ! create or open groups
            name = step // "/parcels"
            step_group = open_h5_group(step)
            group = open_h5_group(name)

            !
            ! write parcel data
            !

            call write_h5_dataset_2d(name, "position", parcels%position(1:n_parcels, :))
            call write_h5_dataset_2d(name, "velocity", parcels%velocity(1:n_parcels, :))

            if (allocated(parcels%stretch)) then
                call write_h5_dataset_1d(name, "stretch", parcels%stretch(1:n_parcels))
            endif

            if (allocated(parcels%B)) then
                call write_h5_dataset_2d(name, "B", parcels%B(1:n_parcels, :))
            endif

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_parcels

        subroutine split(threshold)
            double precision, intent(in) :: threshold


        end subroutine split

        subroutine alloc_parcel_mem(num)
            integer, intent(in) :: num

            allocate(parcels%position(num, 2))
            allocate(parcels%velocity(num, 2))
            allocate(parcels%stretch(num))
            allocate(parcels%B(num, 2))
            allocate(parcels%volume(num))
        end subroutine alloc_parcel_mem

        subroutine dealloc_parcel_mem
            deallocate(parcels%position)
            deallocate(parcels%velocity)

            if (allocated(parcels%stretch)) then
                deallocate(parcels%stretch)
            endif

            if (allocated(parcels%B)) then
                deallocate(parcels%B)
            endif

            deallocate(parcels%volume)
        end subroutine dealloc_parcel_mem

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcel_container
