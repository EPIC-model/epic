module parcel_container
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_1d, &
                       write_h5_dataset_2d, &
                       open_h5_group,       &
                       get_step_group_name
    use ellipse, only : get_angles
    implicit none

    integer :: n_parcels

    type parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            velocity,   &
            stretch,    &
            volume,     &
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        subroutine write_h5_parcels(iter)
            integer, intent(in)           :: iter ! iteration
            integer(hid_t)                :: group
            integer(hid_t)                :: step_group
            character(:), allocatable     :: step
            character(:), allocatable     :: name
            double precision, allocatable :: angle(:)

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
                call write_h5_dataset_1d(name, "stretch", parcels%stretch(1:n_parcels, 1))
            endif

            call write_h5_dataset_1d(name, "volume", parcels%volume(1:n_parcels, 1))

            if (allocated(parcels%B)) then
                call write_h5_dataset_2d(name, "B", parcels%B(1:n_parcels, :))

                angle = get_angles(parcels%B, n_parcels)
                call write_h5_dataset_1d(name, "orientation", angle)
            endif

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_parcels


        ! overwrite parcel n with parcel m
        subroutine parcel_pack(mask)
            logical, intent(in) :: mask(:)

            parcels%position(:, 1) = pack(parcels%position(:, 1), mask)
            parcels%position(:, 2) = pack(parcels%position(:, 2), mask)

            parcels%velocity(:, 1) = pack(parcels%velocity(:, 1), mask)
            parcels%velocity(:, 2) = pack(parcels%velocity(:, 2), mask)

            if (allocated(parcels%stretch)) then
                parcels%stretch(:, 1)  = pack(parcels%stretch(:, 1), mask)
            endif

            parcels%volume(:, 1)  = pack(parcels%volume(:, 1), mask)

            parcels%B(:, 1) = pack(parcels%B(:, 1), mask)
            parcels%B(:, 2) = pack(parcels%B(:, 2), mask)

            ! update parcel number
            n_parcels = n_parcels - count(.not. mask)

        end subroutine parcel_pack


        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(num, 2))
            allocate(parcels%velocity(num, 2))
            allocate(parcels%stretch(num, 1))
            allocate(parcels%B(num, 2))
            allocate(parcels%volume(num, 1))
        end subroutine parcel_alloc


        subroutine parcel_dealloc
            deallocate(parcels%position)
            deallocate(parcels%velocity)

            if (allocated(parcels%stretch)) then
                deallocate(parcels%stretch)
            endif

            if (allocated(parcels%B)) then
                deallocate(parcels%B)
            endif

            deallocate(parcels%volume)
        end subroutine parcel_dealloc

end module parcel_container
