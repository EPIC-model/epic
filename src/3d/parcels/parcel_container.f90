! =============================================================================
! This module stores the parcel data and provides subroutines to write, modify
! allocate and deallocate it.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, extenti, center, lower, upper
    use parcel_ellipsoid, only : parcel_ellipsoid_allocate, parcel_ellipsoid_deallocate
    use mpi_communicator
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    integer :: n_parcels        ! local number of parcels
    integer :: n_total_parcels  ! global number of parcels (over all MPI ranks)

    ! buffer indices to access individual parcel attributes
    integer, protected :: IDX_X_POS,        & ! x-position
                          IDX_Y_POS,        & ! y-position
                          IDX_Z_POS,        & ! z-position
                          IDX_X_VOR,        & ! x-vorticity
                          IDX_Y_VOR,        & ! y-vorticity
                          IDX_Z_VOR,        & ! z-vorticity
                          IDX_B11,          & ! B11 shape matrix element
                          IDX_B12,          & ! B12 shape matrix element
                          IDX_B13,          & ! B13 shape matrix element
                          IDX_B22,          & ! B22 shape matrix element
                          IDX_B23,          & ! B23 shape matrix element
                          IDX_VOL,          & ! volume
                          IDX_BUO,          & ! buoyancy
#ifndef ENABLE_DRY_MODE
                          IDX_HUM,          & ! humidity
#endif
                          IDX_RK4_X_DPOS,   & ! RK4 variable delta x-position
                          IDX_RK4_Y_DPOS,   & ! RK4 variable delta y-position
                          IDX_RK4_Z_DPOX,   & ! RK4 variable delta z-position
                          IDX_RK4_X_DVOR,   & ! RK4 variable delta x-vorticity
                          IDX_RK4_Y_DVOR,   & ! RK4 variable delta y-vorticity
                          IDX_RK4_Z_DVOR,   & ! RK4 variable delta z-vorticity
                          IDX_RK4_DB11,     & ! RK4 variable for B11
                          IDX_RK4_DB12,     & ! RK4 variable for B12
                          IDX_RK4_DB13,     & ! RK4 variable for B13
                          IDX_RK4_DB22,     & ! RK4 variable for B22
                          IDX_RK4_DB23,     & ! RK4 variable for B23
                          IDX_RK4_DUDX,     & ! RK4 variable du/dx
                          IDX_RK4_DUDY,     & ! RK4 variable du/dy
                          IDX_RK4_DVDY,     & ! RK4 variable dv/dy
                          IDX_RK4_DWDX,     & ! RK4 variable dw/dx
                          IDX_RK4_DWDY        ! RK4 variable dw/dy

    ! number of  parcel attributes
    ! (components are counted individually, e.g. position counts as 3 attributes)
    integer, protected :: n_par_attrib

    type parcel_container_type
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            vorticity,  &
            B               ! B matrix entries; ordering:
                            ! B(:, 1) = B11, B(:, 2) = B12, B(:, 3) = B13
                            ! B(:, 4) = B22, B(:, 5) = B23

        double precision, allocatable, dimension(:) :: &
            volume,     &
#ifndef ENABLE_DRY_MODE
            humidity,   &
#endif
            buoyancy


        ! LS-RK4 variables
        double precision, allocatable, dimension(:, :) :: &
            delta_pos,  &
            delta_vor,  &
            strain,     &
            delta_b

    end type parcel_container_type

    type(parcel_container_type) parcels

    private :: set_buffer_indices

    contains

        ! Sets the buffer indices for sending parcels around
        ! (called in parcel_alloc)
        subroutine set_buffer_indices
            integer :: i

            IDX_X_POS = 1   ! x-position
            IDX_Y_POS = 2   ! y-position
            IDX_Z_POS = 3   ! z-position
            IDX_X_VOR = 4   ! x-vorticity
            IDX_Y_VOR = 5   ! y-vorticity
            IDX_Z_VOR = 6   ! z-vorticity
            IDX_B11   = 7   ! B11 shape matrix element
            IDX_B12   = 8   ! B12 shape matrix element
            IDX_B13   = 9   ! B13 shape matrix element
            IDX_B22   = 10  ! B22 shape matrix element
            IDX_B23   = 11  ! B23 shape matrix element
            IDX_VOL   = 12  ! volume
            IDX_BUO   = 13  ! buoyancy

            i = IDX_BUO + 1
#ifndef ENABLE_DRY_MODE
            IDX_HUM  = i
            i = i + 1
#endif

            ! LS-RK4 variables
            IDX_RK4_X_DPOS = i
            IDX_RK4_Y_DPOS = i + 1
            IDX_RK4_Z_DPOX = i + 2
            IDX_RK4_X_DVOR = i + 3
            IDX_RK4_Y_DVOR = i + 4
            IDX_RK4_Z_DVOR = i + 5
            IDX_RK4_DB11 = i + 6
            IDX_RK4_DB12 = i + 7
            IDX_RK4_DB13 = i + 8
            IDX_RK4_DB22 = i + 9
            IDX_RK4_DB23 = i + 10
            IDX_RK4_DUDX = i + 11
            IDX_RK4_DUDY = i + 12
            IDX_RK4_DVDY = i + 13
            IDX_RK4_DWDX = i + 14
            IDX_RK4_DWDY = i + 15

            n_par_attrib = IDX_RK4_DWDY

        end subroutine set_buffer_indices

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two zonal coordinates
        ! across periodic edges
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
        ! @returns delx = x1 - x2
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_delx_across_periodic(x1, x2) result (delx)
#else
        elemental function get_delx_across_periodic(x1, x2) result (delx)
#endif
            double precision, intent(in) :: x1, x2
            double precision             :: delx

            delx = x1 - x2
#ifndef NDEBUG
            if ((x1 < lower(1)) .or. (x2 < lower(1)) .or. &
                (x1 > upper(1)) .or. (x2 > upper(1))) then
                write(*,*) 'point outside domain was fed into get_delx'
                write(*,*) 'x1, x2, lower(1), upper(1)'
                write(*,*) x1, x2, lower(1), upper(1)
            endif
#endif
            ! works across periodic edge
            delx = delx - extent(1) * dble(nint(delx * extenti(1)))
        end function get_delx_across_periodic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two zonal coordinates
        ! @param[in] x1 first zonal position
        ! @param[in] x2 second zonal position
        ! @returns delx = x1 - x2
        elemental function get_delx(x1, x2) result (delx)
            double precision, intent(in) :: x1, x2
            double precision             :: delx
            delx = x1 - x2
        end function get_delx

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two meridional coordinates
        ! across periodic edges
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        ! WARNING input needs to be between lower and upper (see debug statement)
#ifndef NDEBUG
        function get_dely_across_periodic(y1, y2) result (dely)
#else
        elemental function get_dely_across_periodic(y1, y2) result (dely)
#endif
            double precision, intent(in) :: y1, y2
            double precision             :: dely

            dely = y1 - y2
#ifndef NDEBUG
            if ((y1 < lower(2)) .or. (y2 < lower(2)) .or. &
                (y1 > upper(2)) .or. (y2 > upper(2))) then
                write(*,*) 'point outside domain was fed into get_dely'
                write(*,*) 'y1, y2, lower(2), upper(2)'
                write(*,*) y1, y2, lower(2), upper(2)
            endif
#endif
            ! works across periodic edge
            dely = dely - extent(2) * dble(nint(dely * extenti(2)))
        end function get_dely_across_periodic

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the difference between two meridional coordinates
        ! @param[in] y1 first meridional position
        ! @param[in] y2 second meridional position
        ! @returns dely = y1 - y2
        elemental function get_dely(y1, y2) result (dely)
            double precision, intent(in) :: y1, y2
            double precision             :: dely
            dely = y1 - y2
        end function get_dely

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Overwrite parcel n with parcel m
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine parcel_replace(n, m)
            integer, intent(in) :: n, m

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print '(a19, i0, a6, i0)', '    replace parcel ', n, ' with ', m
            endif
#endif

            parcels%position(:, n) = parcels%position(:, m)

            parcels%vorticity(:, n) = parcels%vorticity(:, m)

            parcels%volume(n)  = parcels%volume(m)
            parcels%buoyancy(n) = parcels%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n) = parcels%humidity(m)
#endif
            parcels%B(:, n) = parcels%B(:, m)

        end subroutine parcel_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Allocate parcel memory
        ! @param[in] num number of parcels
        subroutine parcel_alloc(num)
            integer, intent(in) :: num

            call set_buffer_indices

            allocate(parcels%position(3, num))
            allocate(parcels%vorticity(3, num))
            allocate(parcels%B(5, num))
            allocate(parcels%volume(num))
            allocate(parcels%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(parcels%humidity(num))
#endif
            call parcel_ellipsoid_allocate(num)

            ! LS-RK4 variables
            allocate(parcels%delta_pos(3, num))
            allocate(parcels%delta_vor(3, num))
            allocate(parcels%strain(5, num))
            allocate(parcels%delta_b(5, num))
        end subroutine parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate parcel memory
        subroutine parcel_dealloc

            if (.not. allocated(parcels%position)) then
                return
            endif

            n_parcels = 0

            n_total_parcels = 0
            if (comm%size > 1) then
                call mpi_blocking_reduce(n_total_parcels, MPI_SUM)
            endif

            deallocate(parcels%position)
            deallocate(parcels%vorticity)
            deallocate(parcels%B)
            deallocate(parcels%volume)
            deallocate(parcels%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(parcels%humidity)
#endif
            call parcel_ellipsoid_deallocate

            ! LS-RK4 variables
            deallocate(parcels%delta_pos)
            deallocate(parcels%delta_vor)
            deallocate(parcels%strain)
            deallocate(parcels%delta_b)

        end subroutine parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Serialize all parcel attributes into a single buffer
        subroutine parcel_serialize(n, buffer)
            integer,          intent(in)  :: n
            double precision, intent(out) :: buffer(n_par_attrib)

            buffer(IDX_X_POS:IDX_Z_POS) = parcels%position(:, n)
            buffer(IDX_X_VOR:IDX_Z_VOR) = parcels%vorticity(:, n)
            buffer(IDX_B11:IDX_B23)     = parcels%B(:, n)
            buffer(IDX_VOL)             = parcels%volume(n)
            buffer(IDX_BUO)             = parcels%buoyancy(n)
#ifndef ENABLE_DRY_MODE
            buffer(IDX_HUM)             = parcels%humidity(n)
#endif
        end subroutine parcel_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deserialize all parcel attributes from a single buffer
        subroutine parcel_deserialize(n, buffer)
            integer,          intent(in) :: n
            double precision, intent(in) :: buffer(n_par_attrib)

            parcels%position(:, n)  = buffer(IDX_X_POS:IDX_Z_POS)
            parcels%vorticity(:, n) = buffer(IDX_X_VOR:IDX_Z_VOR)
            parcels%B(:, n)         = buffer(IDX_B11:IDX_B23)
            parcels%volume(n)       = buffer(IDX_VOL)
            parcels%buoyancy(n)     = buffer(IDX_BUO)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(n)     = buffer(IDX_HUM)
#endif
        end subroutine parcel_deserialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_pack(pid, num, buffer)
            integer,          intent(in)  :: pid(:)
            integer,          intent(in)  :: num
            double precision, intent(out) :: buffer(:)
            integer                       :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * n_par_attrib
                j = n * n_par_attrib
                call parcel_serialize(pid(n), buffer(i:j))
            enddo
        end subroutine parcel_pack

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_unpack(num, buffer)
            integer,          intent(in) :: num
            double precision, intent(in) :: buffer(:)
            integer                      :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * n_par_attrib
                j = n * n_par_attrib
                call parcel_deserialize(n_parcels + n, buffer(i:j))
            enddo

            n_parcels = n_parcels + num

        end subroutine parcel_unpack

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This algorithm replaces invalid parcels with valid parcels
        ! from the end of the container
        ! @param[in] pid are the parcel indices of the parcels to be deleted
        ! @param[in] n_del is the array size of pid
        ! @pre
        !   - pid must be sorted in ascending order
        !   - pid must be contiguously filled
        !   The above preconditions must be fulfilled so that the
        !   parcel pack algorithm works correctly.
        subroutine parcel_delete(pid, n_del)
            integer, intent(in) :: pid(0:)
            integer, intent(in) :: n_del
            integer             :: k, l, m

            ! l points always to the last valid parcel
            l = n_parcels

            ! k points always to last invalid parcel in pid
            k = n_del

            ! find last parcel which is not invalid
            do while ((k > 0) .and. (l == pid(k)))
                l = l - 1
                k = k - 1
            enddo

            if (l == 0) then
                print *, "Error: All parcels are invalid."
                stop
            endif

            ! replace invalid parcels with the last valid parcel
            m = 1

            do while (m <= k)
                ! invalid parcel; overwrite *pid(m)* with last valid parcel *l*
                call parcel_replace(pid(m), l)

                l = l - 1

                ! find next valid last parcel
                do while ((k > 0) .and. (l == pid(k)))
                    l = l - 1
                    k = k - 1
                enddo

                ! next invalid
                m = m + 1
            enddo

            ! update number of valid parcels
            n_parcels = n_parcels - n_del

        end subroutine parcel_delete

end module parcel_container
