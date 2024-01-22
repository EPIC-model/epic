! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module surface_parcel_container
    use parameters, only : max_num_surf_parcels
    use merge_sort
    use parcel_ops
    implicit none

    integer :: n_top_parcels, n_bot_parcels

    type surface_parcel_container_type
        double precision, allocatable, dimension(:) :: &
            position,   &           ! denotes the start point of a line
            vorticity,  &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy,   &
            volume

        integer, allocatable, dimension(:) :: &
            right                    ! j = right(i) gives the index j of the parcel on the right side

        ! LS-RK-4 arrays
        double precision, allocatable, dimension(:) :: &
            delta_pos, &
            delta_vor      ! vorticity integration


    end type surface_parcel_container_type

    type(surface_parcel_container_type) top_parcels, bot_parcels

    private :: alloc, dealloc

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        function get_surface_parcel_length(n, spar) result(length)
            integer,                             intent(in) :: n
            type(surface_parcel_container_type), intent(in) :: spar
            integer                                         :: j
            double precision                                :: length

            j = spar%right(n)

            length = get_delx(spar%position(j), spar%position(n))

        end function get_surface_parcel_length

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine surface_parcel_sort(n_par, sp)
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: sp
            integer                                            :: indx(n_par)
            type(surface_parcel_container_type)                :: tmp
            integer                                            :: n, i

            call alloc(max_num_surf_parcels, tmp)

            tmp%position = sp%position!(1:n_par)

            ! sort position in ascending order
            call msort(tmp%position(1:n_par), indx)

            do n = 1, n_par
                i = indx(n)

                tmp%vorticity(n) = sp%vorticity(i)
                tmp%buoyancy(n) = sp%buoyancy(i)
                tmp%volume(n) = sp%volume(i)
#ifndef ENABLE_DRY_MODE
                tmp%humidity(n) = sp%humidity(i)
#endif

                tmp%right(n) = mod(n, n_par) + 1

            enddo

            call move_alloc(from=tmp%position, to=sp%position)
            call move_alloc(from=tmp%vorticity, to=sp%vorticity)
            call move_alloc(from=tmp%buoyancy, to=sp%buoyancy)
            call move_alloc(from=tmp%volume, to=sp%volume)
#ifndef ENABLE_DRY_MODE
            call move_alloc(from=tmp%humidity, to=sp%humidity)
#endif
            call move_alloc(from=tmp%right, to=sp%right)

        end subroutine surface_parcel_sort

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine surface_parcel_replace(n, m, sp)
            integer,                             intent(in)    :: n, m
            type(surface_parcel_container_type), intent(inout) :: sp

            sp%position(n) = sp%position(m)

            sp%vorticity(n) = sp%vorticity(m)

            sp%buoyancy(n) = sp%buoyancy(m)
            sp%volume(n) = sp%volume(m)
#ifndef ENABLE_DRY_MODE
            sp%humidity(n) = sp%humidity(m)
#endif
            sp%right(n) = sp%right(m)
        end subroutine surface_parcel_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine alloc(num, sp)
            integer,                             intent(in)    :: num
            type(surface_parcel_container_type), intent(inout) :: sp

            allocate(sp%position(num))
            allocate(sp%vorticity(num))
            allocate(sp%buoyancy(num))
            allocate(sp%volume(num))
#ifndef ENABLE_DRY_MODE
            allocate(sp%humidity(num))
#endif
            allocate(sp%right(num))

            allocate(sp%delta_pos(num))
            allocate(sp%delta_vor(num))

        end subroutine alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine dealloc(sp)
            type(surface_parcel_container_type), intent(inout) :: sp

            deallocate(sp%position)
            deallocate(sp%vorticity)
            deallocate(sp%buoyancy)
            deallocate(sp%volume)
#ifndef ENABLE_DRY_MODE
            deallocate(sp%humidity)
#endif
            deallocate(sp%right)

            deallocate(sp%delta_pos)
            deallocate(sp%delta_vor)
        end subroutine dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine surface_parcel_alloc(num)
            integer, intent(in) :: num

            call alloc(num, top_parcels)
            call alloc(num, bot_parcels)

        end subroutine surface_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate parcel memory
        subroutine surface_parcel_dealloc

            call dealloc(top_parcels)
            call dealloc(bot_parcels)

        end subroutine surface_parcel_dealloc

end module surface_parcel_container
