! =============================================================================
! This module provides the base class for parcel containers.
! =============================================================================
module parcel_container
    use options, only : verbose
    use parameters, only : extent, extenti, center, lower, upper
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error
    use armanip, only : resize_array
    use mpi_timer, only : start_timer, stop_timer
    implicit none

    integer :: resize_timer

    ! Parcel container type
    type, abstract :: pc_type

            ! number of  parcel attributes
            ! (components are counted individually, e.g. position counts as 3 attributes)
            integer :: attr_num     ! number of parcel attributes
            integer :: mix_attr_num ! number of parcel attributes used in mixing
            integer :: local_num    ! local number of parcels
            integer :: total_num    ! global number of parcels (over all MPI ranks)
            integer :: max_num      ! capacity per attribute, i.e. maximum number of parcels
            integer :: nz           ! number of grid cells this parcel container contributes to
                                    ! (used for the MPI communication)

            ! ---------------------------------------------------------------------
            !   Parcel attributes (common to all types):
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

            ! ---------------------------------------------------------------------
            ! LS-RK variables:
            double precision, allocatable, dimension(:, :) :: &
                delta_pos,  &   ! velocity
                delta_vor,  &   ! vorticity tendency
                strain,     &
                delta_b         ! B-matrix tendency

            ! -------------------------------------------------------------------------
            ! buffer indices to access parcel attributes for (de-)serialization;
            ! the buffer indices are set in the extendend types
            integer :: IDX_POS_BEG,         & ! position vector begin
                       IDX_POS_END,         & ! position vector end
                       IDX_VOR_BEG,         & ! vorticity vector begin
                       IDX_VOR_END,         & ! vorticity vector end
                       IDX_SHAPE_BEG,       & ! shape matrix begin
                       IDX_SHAPE_END,       & ! shape matrix end
                       IDX_VOL,             & ! volume
                       IDX_BUO,             & ! buoyancy
#ifndef ENABLE_DRY_MODE
                       IDX_HUM,             & ! humidity
#endif
                       IDX_RK_POS_BEG,      & ! RK variable delta position vector begin
                       IDX_RK_POS_END,      & ! RK variable delta position vecto end
                       IDX_RK_VOR_BEG,      & ! RK variable delta vorticity vector begin
                       IDX_RK_VOR_END,      & ! RK variable delta vorticity vector end
                       IDX_RK_SHAPE_BEG,    & ! RK variable shape vector begin
                       IDX_RK_SHAPE_END,    & ! RK variable shape vector end
                       IDX_RK_STRAIN_BEG,   & ! RK variable strain vector begin
                       IDX_RK_STRAIN_END      ! RK variable strain vector end

        contains
            ! Base procedures (usually called in derived class procedures):
            procedure :: parcel_base_allocate
            procedure :: parcel_base_deallocate
            procedure :: parcel_base_replace
            procedure :: parcel_base_resize
            procedure :: parcel_base_serialize
            procedure :: parcel_base_deserialize
            ! Basic procedures common to all derived classes:
            procedure :: mixing_serialize     => parcel_mixing_serialize
            procedure :: mixing_deserialize   => parcel_mixing_deserialize
            procedure :: pack                 => parcel_pack
            procedure :: unpack               => parcel_unpack
            procedure :: delete               => parcel_delete
            ! Pure virtual procedures:
            procedure(parcel_allocate),       deferred :: allocate
            procedure(parcel_deallocate),     deferred :: deallocate
            procedure(parcel_serialize),      deferred :: serialize
            procedure(parcel_deserialize),    deferred :: deserialize
            procedure(parcel_replace),        deferred :: replace
            procedure(parcel_resize),         deferred :: resize
            procedure(parcel_is_small),       deferred :: is_small
            procedure(parcel_get_z_position), deferred :: get_z_position
            procedure(parcel_get_position),   deferred :: get_position
            procedure(parcel_set_position),   deferred :: set_position
            procedure(parcel_set_volume),     deferred :: set_volume

    end type pc_type

        interface
            subroutine parcel_allocate(this, num)
                import :: pc_type
                class(pc_type), intent(inout) :: this
                integer,        intent(in)    :: num
            end subroutine parcel_allocate

            subroutine parcel_deallocate(this)
                import :: pc_type
                class(pc_type), intent(inout) :: this
            end subroutine parcel_deallocate

            subroutine parcel_replace(this, n, m)
                import :: pc_type
                class(pc_type), intent(inout) :: this
                integer,        intent(in)    :: n, m
            end subroutine

            subroutine parcel_resize(this, new_size)
                import :: pc_type
                class(pc_type), intent(inout) :: this
                integer,        intent(in)    :: new_size
            end subroutine parcel_resize

            subroutine parcel_serialize(this, n, buffer)
                import :: pc_type
                class(pc_type),   intent(in)  :: this
                integer,          intent(in)  :: n
                double precision, intent(out) :: buffer(this%attr_num)
            end subroutine parcel_serialize

            subroutine parcel_deserialize(this, n, buffer)
                import :: pc_type
                class(pc_type),   intent(inout) :: this
                integer,          intent(in)    :: n
                double precision, intent(in)    :: buffer(this%attr_num)
            end subroutine parcel_deserialize

            logical pure function parcel_is_small(this, n)
                import :: pc_type
                class(pc_type), intent(in) :: this
                integer,        intent(in) :: n
            end function parcel_is_small

            double precision pure function parcel_get_z_position(this, n)
                import :: pc_type
                class(pc_type),   intent(in)  :: this
                integer,          intent(in)  :: n
            end function parcel_get_z_position

            pure function parcel_get_position(this, n) result(pos)
                import :: pc_type
                class(pc_type), intent(in) :: this
                integer,        intent(in) :: n
                double precision           :: pos(3)
            end function parcel_get_position

            subroutine parcel_set_position(this, n, pos)
                import :: pc_type
                class(pc_type),   intent(inout) :: this
                integer,          intent(in)    :: n
                double precision, intent(in)    :: pos(3)
            end subroutine parcel_set_position

            subroutine parcel_set_volume(this, n, vol)
                import :: pc_type
                class(pc_type),   intent(inout) :: this
                integer,          intent(in)    :: n
                double precision, intent(in)    :: vol
            end subroutine parcel_set_volume
        end interface

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Allocate parcel memory
        ! ATTENTION: Extended types must allocate additional parcel attributes
        !            in their own routine.
        ! @param[in] num number of parcels
        ! @param[in] n_pos number of spatial dimensions
        ! @param[in] n_vec number of dimensions of vector attributes
        ! @param[in] n_shape number of B matrix elements
        ! @param[in] n_strain number of strain elements
        subroutine parcel_base_allocate(this, num, n_pos, n_vec, n_shape, n_strain)
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: num
            integer,        intent(in)    :: n_pos, n_vec, n_shape, n_strain

            allocate(this%position(n_pos, num))
            allocate(this%vorticity(n_pos, num))
            allocate(this%B(n_shape, num))
            allocate(this%volume(num))
            allocate(this%buoyancy(num))
#ifndef ENABLE_DRY_MODE
            allocate(this%humidity(num))
#endif
            ! LS-RK variables
            allocate(this%delta_pos(n_pos, num))
            allocate(this%delta_vor(n_pos, num))
            allocate(this%strain(n_strain, num))
            allocate(this%delta_b(n_shape, num))

            this%max_num = num

            ! vorticity, buoyancy, volume
            this%mix_attr_num = n_pos + 2
#ifndef ENABLE_DRY_MODE
            ! and humidity
            this%mix_attr_num = this%mix_attr_num + 1
#endif

        end subroutine parcel_base_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate parcel memory
        ! ATTENTION: Extended types must deallocate additional parcel attributes
        !            in their own routine.
        subroutine parcel_base_deallocate(this)
            class(pc_type), intent(inout) :: this

            if (.not. allocated(this%position)) then
                return
            endif

            this%local_num = 0
            this%total_num = 0
            this%max_num   = 0

            deallocate(this%position)
            deallocate(this%vorticity)
            deallocate(this%B)
            deallocate(this%volume)
            deallocate(this%buoyancy)
#ifndef ENABLE_DRY_MODE
            deallocate(this%humidity)
#endif
            ! LS-RK variables
            deallocate(this%delta_pos)
            deallocate(this%delta_vor)
            deallocate(this%strain)
            deallocate(this%delta_b)

        end subroutine parcel_base_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Overwrite parcel n with parcel m. This subroutine only replaces the
        ! common types.
        ! ATTENTION: Extended types must replace additional parcel attributes
        !            in their own routine.
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine parcel_base_replace(this, n, m)
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: n, m

            this%position(:, n)  = this%position(:, m)
            this%vorticity(:, n) = this%vorticity(:, m)
            this%volume(n)       = this%volume(m)
            this%buoyancy(n)     = this%buoyancy(m)
#ifndef ENABLE_DRY_MODE
            this%humidity(n)     = this%humidity(m)
#endif
            this%B(:, n)         = this%B(:, m)

            ! LS-RK variables:
            this%delta_pos(:, n) = this%delta_pos(:, m)
            this%delta_vor(:, n) = this%delta_vor(:, m)
            this%delta_b(:, n)   = this%delta_b(:, m)
            this%strain(:, n)    = this%strain(:, m)

        end subroutine parcel_base_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Resize the parcel container
        ! ATTENTION: Extended types must resize additional parcel attributes
        !            in their own routine.
        ! @param[in] new_size is the new size of each attribute
        subroutine parcel_base_resize(this, new_size)
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call start_timer(resize_timer)

            if (new_size < this%local_num) then
                call mpi_exit_on_error(&
                    "in parcel_container::parcel_base_resize: losing parcels when resizing.")
            endif

            this%max_num = new_size

            call resize_array(this%position, new_size, this%local_num)

            call resize_array(this%vorticity, new_size, this%local_num)
            call resize_array(this%B, new_size, this%local_num)
            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%buoyancy, new_size, this%local_num)
#ifndef ENABLE_DRY_MODE
            call resize_array(this%humidity, new_size, this%local_num)
#endif

            ! LS-RK variables
            call resize_array(this%delta_pos, new_size, this%local_num)
            call resize_array(this%delta_vor, new_size, this%local_num)
            call resize_array(this%strain, new_size, this%local_num)
            call resize_array(this%delta_b, new_size, this%local_num)

            call stop_timer(resize_timer)

        end subroutine parcel_base_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Serialize all parcel attributes into a single buffer
        subroutine parcel_base_serialize(this, n, buffer)
            class(pc_type),   intent(in)  :: this
            integer,          intent(in)  :: n
            double precision, intent(out) :: buffer(this%attr_num)

            buffer(this%IDX_POS_BEG:this%IDX_POS_END)       = this%position(:, n)
            buffer(this%IDX_VOR_BEG:this%IDX_VOR_END)       = this%vorticity(:, n)
            buffer(this%IDX_SHAPE_BEG:this%IDX_SHAPE_END)   = this%B(:, n)
            buffer(this%IDX_VOL)                            = this%volume(n)
            buffer(this%IDX_BUO)                            = this%buoyancy(n)
#ifndef ENABLE_DRY_MODE
            buffer(this%IDX_HUM)                            = this%humidity(n)
#endif
            ! LS-RK variables:
            buffer(this%IDX_RK_POS_BEG:this%IDX_RK_POS_END)       = this%delta_pos(:, n)
            buffer(this%IDX_RK_VOR_BEG:this%IDX_RK_VOR_END)       = this%delta_vor(:, n)
            buffer(this%IDX_RK_SHAPE_BEG:this%IDX_RK_SHAPE_END)   = this%delta_b(:, n)
            buffer(this%IDX_RK_STRAIN_BEG:this%IDX_RK_STRAIN_END) = this%strain(:, n)

        end subroutine parcel_base_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deserialize all parcel attributes from a single buffer
        subroutine parcel_base_deserialize(this, n, buffer)
            class(pc_type),   intent(inout) :: this
            integer,          intent(in)    :: n
            double precision, intent(in)    :: buffer(this%attr_num)

            this%position(:, n)  = buffer(this%IDX_POS_BEG:this%IDX_POS_END)
            this%vorticity(:, n) = buffer(this%IDX_VOR_BEG:this%IDX_VOR_END)
            this%B(:, n)         = buffer(this%IDX_SHAPE_BEG:this%IDX_SHAPE_END)
            this%volume(n)       = buffer(this%IDX_VOL)
            this%buoyancy(n)     = buffer(this%IDX_BUO)
#ifndef ENABLE_DRY_MODE
            this%humidity(n)     = buffer(this%IDX_HUM)
#endif
            ! LS-RK variables:
            this%delta_pos(:, n) = buffer(this%IDX_RK_POS_BEG:this%IDX_RK_POS_END)
            this%delta_vor(:, n) = buffer(this%IDX_RK_VOR_BEG:this%IDX_RK_VOR_END)
            this%delta_b(:, n)   = buffer(this%IDX_RK_SHAPE_BEG:this%IDX_RK_SHAPE_END)
            this%strain(:, n)    = buffer(this%IDX_RK_STRAIN_BEG:this%IDX_RK_STRAIN_END)

        end subroutine parcel_base_deserialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Serialize all parcel attributes to mix into a single buffer
        subroutine parcel_mixing_serialize(this, n, buffer)
            class(pc_type),   intent(in)  :: this
            integer,          intent(in)  :: n
            double precision, intent(out) :: buffer(this%mix_attr_num)

            buffer(1:2) = this%vorticity(1:2, n)
            buffer(3)   = this%volume(n)
            buffer(4)   = this%buoyancy(n)
#ifndef ENABLE_DRY_MODE
            buffer(5)   = this%humidity(n)
#endif
        end subroutine parcel_mixing_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deserialize all parcel attributes to mix from a single buffer
        subroutine parcel_mixing_deserialize(this, n, buffer)
            class(pc_type),   intent(inout) :: this
            integer,          intent(in)    :: n
            double precision, intent(in)    :: buffer(this%mix_attr_num)

            this%vorticity(1:2, n) = buffer(1:2)
            this%volume(n)       = buffer(3)
            this%buoyancy(n)     = buffer(4)
#ifndef ENABLE_DRY_MODE
            this%humidity(n)     = buffer(5)
#endif
        end subroutine parcel_mixing_deserialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_pack(this, pid, num, buffer)
            class(pc_type),   intent(in)  :: this
            integer,          intent(in)  :: pid(:)
            integer,          intent(in)  :: num
            double precision, intent(out) :: buffer(:)
            integer                       :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * this%attr_num
                j = n * this%attr_num
                call this%serialize(pid(n), buffer(i:j))
            enddo
        end subroutine parcel_pack

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_unpack(this, num, buffer)
            class(pc_type),   intent(inout) :: this
            integer,          intent(in)    :: num
            double precision, intent(in)    :: buffer(:)
            integer                         :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * this%attr_num
                j = n * this%attr_num
                call this%deserialize(this%local_num + n, buffer(i:j))
            enddo

            this%local_num = this%local_num + num

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
        subroutine parcel_delete(this, pid, n_del)
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: pid(0:)
            integer,        intent(in)    :: n_del
            integer                       :: k, l, m

            ! l points always to the last valid parcel
            l = this%local_num

            ! k points always to last invalid parcel in pid
            k = n_del

            ! find last parcel which is not invalid
            do while ((k > 0) .and. (l == pid(k)))
                l = l - 1
                k = k - 1
            enddo

            if (l == -1) then
                call mpi_exit_on_error(&
                    "in parcel_container::parcel_delete: more than all parcels are invalid.")
            endif

            ! replace invalid parcels with the last valid parcel
            m = 1

            do while (m <= k)
                ! invalid parcel; overwrite *pid(m)* with last valid parcel *l*
                call this%replace(pid(m), l)

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
            this%local_num = this%local_num - n_del

        end subroutine parcel_delete

end module parcel_container
