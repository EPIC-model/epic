! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use hdf5
    use writer, only : h5file,              &
                       h5err,               &
                       write_h5_dataset_1d, &
                       write_h5_dataset_2d, &
                       open_h5_group,       &
                       get_step_group_name
    use options, only : verbose
    use parameters, only : extent, hli
    use ellipse, only : get_angles
    implicit none

    integer :: n_parcels

    type parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            velocity,   &
            stretch,    &
            B               ! B matrix entries; ordering B(:, 1) = B11, B(:, 2) = B12

        double precision, allocatable, dimension(:) :: &
            volume,     &
            vorticity,  &
            buoyancy,   &
            humidity
    end type parcel_container_type

    type(parcel_container_type) parcels


    contains

        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
            ! works across periodic edge
            delx = delx - extent(1) * dble(int(delx * hli(1)))
        end function get_delx

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
            call write_h5_dataset_1d(name, "vorticity", parcels%vorticity(1:n_parcels))

            if (allocated(parcels%stretch)) then
                call write_h5_dataset_1d(name, "stretch", parcels%stretch(1:n_parcels, 1))
            endif

            call write_h5_dataset_1d(name, "volume", parcels%volume(1:n_parcels))

            call write_h5_dataset_1d(name, "buoyancy", parcels%buoyancy(1:n_parcels))

            call write_h5_dataset_1d(name, "humidity", parcels%humidity(1:n_parcels))

            if (allocated(parcels%B)) then
                call write_h5_dataset_2d(name, "B", parcels%B(1:n_parcels, :))

                angle = get_angles(parcels%B, parcels%volume, n_parcels)
                call write_h5_dataset_1d(name, "orientation", angle)
            endif

            call h5gclose_f(group, h5err)
            call h5gclose_f(step_group, h5err)

        end subroutine write_h5_parcels


        ! overwrite parcel n with parcel m
        subroutine parcel_replace(n, m)
            integer, intent(in) :: n, m

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print '(a19, i0, a6, i0)', '    replace parcel ', n, ' with ', m
            endif
#endif

            parcels%position(n, 1) = parcels%position(m, 1)
            parcels%position(n, 2) = parcels%position(m, 2)

            parcels%velocity(n, 1) = parcels%velocity(m, 1)
            parcels%velocity(n, 2) = parcels%velocity(m, 2)

            parcels%vorticity(n) = parcels%vorticity(m)

            if (allocated(parcels%stretch)) then
                parcels%stretch(n, 1)  = parcels%stretch(m, 1)
            endif

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
            parcels%humidity(n) = parcels%humidity(m)

            parcels%B(n, 1) = parcels%B(m, 1)
            parcels%B(n, 2) = parcels%B(m, 2)

        end subroutine parcel_replace


        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            allocate(parcels%position(num, 2))
            allocate(parcels%velocity(num, 2))
            allocate(parcels%vorticity(num))
            allocate(parcels%stretch(num, 1))
            allocate(parcels%B(num, 2))
            allocate(parcels%volume(num))
            allocate(parcels%buoyancy(num))
            allocate(parcels%humidity(num))
        end subroutine parcel_alloc


        subroutine parcel_dealloc
            deallocate(parcels%position)
            deallocate(parcels%velocity)
            deallocate(parcels%vorticity)

            if (allocated(parcels%stretch)) then
                deallocate(parcels%stretch)
            endif

            if (allocated(parcels%B)) then
                deallocate(parcels%B)
            endif

            deallocate(parcels%volume)
            deallocate(parcels%buoyancy)
            deallocate(parcels%humidity)
        end subroutine parcel_dealloc

end module parcel_container
