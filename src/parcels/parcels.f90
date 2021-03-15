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
            B11, B12,   & ! B matrix entries
            volume

        double precision, allocatable, dimension(:, :) :: &
            position,   &
            velocity
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

            if (allocated(parcels%B11) .and. allocated(parcels%B12)) then
                call write_h5_dataset_1d(name, "B11", parcels%B11(1:n_parcels))
                call write_h5_dataset_1d(name, "B12", parcels%B12(1:n_parcels))
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
            allocate(parcels%B11(num))
            allocate(parcels%B12(num))
            allocate(parcels%volume(num))
        end subroutine alloc_parcel_mem

        subroutine dealloc_parcel_mem
            deallocate(parcels%position)
            deallocate(parcels%velocity)

            if (allocated(parcels%stretch)) then
                deallocate(parcels%stretch)
            endif

            if (allocated(parcels%B11)) then
                deallocate(parcels%B11)
                deallocate(parcels%B12)
            endif

            deallocate(parcels%volume)
        end subroutine dealloc_parcel_mem

        subroutine create(num)
            integer, intent(in) :: num

        end subroutine create

end module parcel_container
