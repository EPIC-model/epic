module parcel_container
    use datatypes, only : int64
    use armanip, only : resize_array
    implicit none

    ! TODO for EPIC implementation
    ! Replace commands by MPI versions

    integer :: resize_timer

    type attr_ptr
        double precision, pointer :: aptr(:)
        character(len=32) :: name
        character(len=128) :: long_name
        character(len=128) :: std_name
        character(len=32) :: dtype ! will need to be MPI_Datatype
        character(len=32) :: unit
        logical :: write_to_netcdf = .false.
    end type attr_ptr

    type int_attr_ptr
        integer(kind=8), pointer :: aptr(:)
        character(len=32) :: name
        character(len=128) :: long_name
        character(len=128) :: std_name
        character(len=32) :: dtype ! will need to be MPI_Datatype
        character(len=32) :: unit
        logical :: write_to_netcdf = .false.
    end type int_attr_ptr

    type, abstract :: base_parcel_type
        double precision, allocatable, dimension(:,:) :: position
        double precision, allocatable, dimension(:,:) :: delta_pos ! For time-integration of the position
                                                                   ! Sedimentation may need to be done separately
        type(attr_ptr), allocatable, dimension(:) :: attrib
        type(int_attr_ptr), allocatable, dimension(:) :: int_attrib
        integer             :: attr_num     ! number of parcel attributes for serialisation
        integer             :: int_attr_num ! number of integer(kind=8) parcel attributes for serialisation
        integer             :: local_num    ! local number of parcels
        integer(kind=int64) :: total_num    ! global number of parcels (over all MPI ranks)
        integer             :: max_num      ! capacity per attribute, i.e. maximum number of parcels
        integer, private :: n_pos = -1      ! number of spatial dimensions

        contains
            procedure :: base_alloc   => base_parcel_alloc
            procedure :: base_dealloc => base_parcel_dealloc
            procedure :: base_resize  => base_parcel_resize
            procedure :: serialize    => base_parcel_serialize
            procedure :: deserialize  => base_parcel_deserialize
            procedure :: replace      => base_parcel_replace
            procedure :: pack         => parcel_pack
            procedure :: unpack       => parcel_unpack
            procedure :: delete       => parcel_delete
            procedure(parcel_alloc), deferred :: alloc
            procedure(parcel_dealloc), deferred :: dealloc
            procedure(parcel_resize), deferred :: resize
            procedure :: print_me
            procedure :: set_dimension
            procedure :: register_attribute
            procedure :: register_int_attribute
            procedure :: reset_attribute
            procedure :: reset_int_attribute

    end type

    ! This type is for parcels associated with dynamics (which will always have volume and voriticity in EPIC)
    type, extends(base_parcel_type) :: dynamic_parcel_type ! add procedures
        double precision, allocatable, dimension(:)   :: volume
        double precision, allocatable, dimension(:,:) :: vorticity
        double precision, allocatable, dimension(:,:) :: delta_vor
        double precision, allocatable, dimension(:) :: dilution
        integer(kind=8), allocatable, dimension(:) :: label
        logical   :: has_labels = .false.

        contains
            procedure :: alloc => dynamic_parcel_alloc
            procedure :: dealloc => dynamic_parcel_dealloc
            procedure :: resize => dynamic_parcel_resize

    end type

    interface
        subroutine parcel_alloc(this, num)
            import base_parcel_type
            class(base_parcel_type), intent(inout) :: this
            integer,                 intent(in)    :: num
        end subroutine parcel_alloc

        subroutine parcel_dealloc(this)
            import base_parcel_type
            class(base_parcel_type), intent(inout) :: this
        end subroutine parcel_dealloc

        subroutine parcel_resize(this, new_size)
            import base_parcel_type
            class(base_parcel_type), intent(inout) :: this
            integer,                 intent(in)    :: new_size
        end subroutine parcel_resize

        logical pure function parcel_is_small(this, n)
            import :: base_parcel_type
            class(base_parcel_type), intent(in) :: this
            integer,        intent(in) :: n
        end function parcel_is_small

    end interface

    interface try_deallocate
        module procedure :: try_deallocate_1d
        module procedure :: try_deallocate_1d_integer
        module procedure :: try_deallocate_2d
    end interface try_deallocate

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine try_deallocate_1d(in_array)
            double precision, allocatable, dimension(:) :: in_array

            if (allocated(in_array)) then
                deallocate(in_array)
            endif

        end subroutine try_deallocate_1d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine try_deallocate_1d_integer(in_array)
            integer(kind=8), allocatable, dimension(:) :: in_array

            if (allocated(in_array)) then
                deallocate(in_array)
            endif

        end subroutine try_deallocate_1d_integer

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine try_deallocate_2d(in_array)
            double precision, allocatable, dimension(:,:) :: in_array

            if (allocated(in_array)) then
                deallocate(in_array)
            endif

        end subroutine try_deallocate_2d

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Allocate parcel memory
        ! ATTENTION: Extended types must allocate additional parcel attributes
        !            in their own routine.
        ! @param[in] num number of parcels
        subroutine base_parcel_alloc(this, num)
            class(base_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            character(len=1)              :: dir(3)
            integer                       :: n

            this%max_num = num
            this%local_num = num ! Usually this is done in parcel_init
            this%attr_num = 0
            this%int_attr_num = 0

            if (this%n_pos > 3) then
                print *, "Only 3 dimensions allowed."
                stop
            endif

            dir = (/'x', 'y', 'z'/)

            allocate(this%position(this%n_pos, num))
            allocate(this%delta_pos(this%n_pos, num))

            do n = 1, this%n_pos
                call this%register_attribute(this%position(n, :), dir(n) // "_position", "m")
                call this%register_attribute(this%delta_pos(n, :), dir(n) // "_position_rk_tendency", "m/s")
            enddo

        end subroutine base_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deallocate parcel memory
        ! ATTENTION: Extended types must deallocate additional parcel attributes
        !            in their own routine.
        subroutine base_parcel_dealloc(this)
            class(base_parcel_type), intent(inout) :: this

            this%local_num = 0
            this%total_num = 0
            this%max_num   = 0
            this%attr_num   = 0
            this%int_attr_num   = 0

            call try_deallocate(this%position)
            call try_deallocate(this%delta_pos)

            if (allocated(this%attrib)) then
                deallocate(this%attrib)
            endif

            if (allocated(this%int_attrib)) then
                deallocate(this%int_attrib)
            endif

        end subroutine base_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine base_parcel_resize(this, new_size)
            class(base_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: new_size
            character(len=1)              :: dir(3)
            integer                       :: n

            if (new_size < this%local_num) then
                print *, "in parcel_container::base_parcel_resize: losing parcels when resizing."
                stop
            endif

            this%max_num = new_size

            if (this%n_pos > 3) then
                print *, "Only 3 dimensions allowed."
                stop
            endif

            dir = (/'x', 'y', 'z'/)

            call resize_array(this%position, new_size, this%local_num)
            call resize_array(this%delta_pos, new_size, this%local_num)

            do n = 1, this%n_pos
                call this%reset_attribute(this%position(n, :), dir(n) // "_position")
                call this%reset_attribute(this%delta_pos(n, :), dir(n) // "_position_rk_tendency")
            enddo

        end subroutine base_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_dimension(this, n_dim)
            class(base_parcel_type), intent(inout) :: this
            integer,                 intent(in)    :: n_dim

            this%n_pos = n_dim

            if (this%n_pos /= -1) then
                print *, "WARNING: Dimension already set."
                return
            endif

            if ((n_dim < 1) .or. (n_dim > 3)) then
                print *, "Only 1-, 2- or 3-dimensional."
                stop
            endif

        end subroutine set_dimension

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine register_attribute(this, attr, name, unit)
            class(base_parcel_type), intent(inout) :: this
            double precision, dimension(:), target, intent(inout) :: attr
            character(len=*) :: name
            character(len=*) :: unit
            type(attr_ptr), allocatable, dimension(:) :: tmp
            integer :: n_attr, n, j

            attr=0.0 ! set to zero when registring, for safety

            ! check if name is unique
            do j = 1, this%attr_num
                 if(trim(this%attrib(j)%name) == trim(name)) then
                     print *, "Attribute name not unique"
                     stop
                 end if
            end do

            do j = 1, this%int_attr_num
                 if(trim(this%int_attrib(j)%name) == trim(name)) then
                     print *, "Attribute name not unique"
                     stop
                 end if
            end do

            this%attr_num = this%attr_num + 1

            if (.not. allocated(this%attrib)) then
                allocate(this%attrib(1))
                this%attrib(1)%aptr => attr
                this%attrib(1)%name = name
                this%attrib(1)%unit = unit
            else
                n_attr = size(this%attrib)

                allocate(tmp(n_attr+1))

                do n = 1, n_attr
                    tmp(n)%aptr => this%attrib(n)%aptr
                    tmp(n)%name = this%attrib(n)%name
                    tmp(n)%unit = this%attrib(n)%unit
                enddo

                tmp(n_attr+1)%aptr => attr
                tmp(n_attr+1)%name = name
                tmp(n_attr+1)%unit = unit

                deallocate(this%attrib)

                call move_alloc(from=tmp, to=this%attrib)

            endif

        end subroutine register_attribute

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reset_attribute(this, attr, name)
            class(base_parcel_type), intent(inout) :: this
            double precision, dimension(:), target, intent(inout) :: attr
            character(len=*) :: name
            integer :: j
            logical :: l_success

            l_success= .false.

             ! check if name is unique
            do j = 1, this%attr_num
                 if(trim(this%attrib(j)%name) == trim(name)) then
                    this%attrib(j)%aptr => attr
                    l_success=.true.
                    exit
                 end if
            end do

            if(.not. l_success) then
                print *, "Attribute not found."
                stop
            end if

        end subroutine reset_attribute

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine register_int_attribute(this, int_attr, name, unit)
            class(base_parcel_type), intent(inout) :: this
            integer(kind=8), dimension(:), target, intent(inout) :: int_attr
            character(len=*) :: name
            character(len=*) :: unit
            type(int_attr_ptr), allocatable, dimension(:) :: tmp
            integer :: n_int_attr, n, j

            int_attr=0 ! set to zero when registring, for safety

            ! check if name is unique
            do j = 1, this%attr_num
                 if(trim(this%attrib(j)%name) == trim(name)) then
                     print *, "Attribute name not unique"
                     stop
                 end if
            end do

            do j = 1, this%int_attr_num
                 if(trim(this%int_attrib(j)%name) == trim(name)) then
                     print *, "Attribute name not unique"
                     stop
                 end if
            end do

            this%int_attr_num = this%int_attr_num + 1

            n_int_attr = size(this%int_attrib)

            if (.not. allocated(this%int_attrib)) then
                allocate(this%int_attrib(1))
                this%int_attrib(1)%aptr => int_attr
                this%int_attrib(1)%name = name
                this%int_attrib(1)%unit = unit
            else
                allocate(tmp(n_int_attr+1))

                do n = 1, n_int_attr
                    tmp(n)%aptr => this%int_attrib(n)%aptr
                    tmp(n)%name = this%int_attrib(n)%name
                    tmp(n)%unit = this%int_attrib(n)%unit
                enddo

                tmp(n_int_attr+1)%aptr => int_attr
                tmp(n_int_attr+1)%name = name
                tmp(n_int_attr+1)%unit = unit

                deallocate(this%int_attrib)

                call move_alloc(from=tmp, to=this%int_attrib)

            endif

        end subroutine register_int_attribute

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine reset_int_attribute(this, int_attr, name)
            class(base_parcel_type), intent(inout) :: this
            integer(kind=8), dimension(:), target, intent(inout) :: int_attr
            character(len=*) :: name
            integer :: j
            logical :: l_success

            l_success= .false.

             ! check if name is unique
            do j = 1, this%int_attr_num
                 if(trim(this%int_attrib(j)%name) == trim(name)) then
                     this%int_attrib(j)%aptr => int_attr
                     l_success=.true.
                     exit
                 end if
            end do

            if(.not. l_success) then
                print *, "Attribute not found."
                stop
            end if

        end subroutine reset_int_attribute

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine dynamic_parcel_alloc(this, num)
            class(dynamic_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%base_alloc(num)

            allocate(this%volume(num))
            call this%register_attribute(this%volume, "volume", "m^3")

            if ((this%n_pos < 2) .or. (this%n_pos > 3)) then
                print *, "Only 2- or 3-dimensional."
                stop
            endif

            if (this%n_pos == 2) then
                allocate(this%vorticity(1, num))
                call this%register_attribute(this%vorticity(1, :), "z_vorticity", "1/s")
                allocate(this%delta_vor(1, num))
                call this%register_attribute(this%delta_vor(1, :), "z_vorticity_rk_tendency", "1/s")
            elseif  (this%n_pos == 3) then
                allocate(this%vorticity(3, num))
                call this%register_attribute(this%vorticity(1, :), "x_vorticity", "1/s")
                call this%register_attribute(this%vorticity(2, :), "y_vorticity", "1/s")
                call this%register_attribute(this%vorticity(3, :), "z_vorticity", "1/s")
                allocate(this%delta_vor(3, num))
                call this%register_attribute(this%delta_vor(1, :), "x_vorticity_rk_tendency", "1/s")
                call this%register_attribute(this%delta_vor(2, :), "y_vorticity_rk_tendency", "1/s")
                call this%register_attribute(this%delta_vor(3, :), "z_vorticity_rk_tendency", "1/s")
            endif

            if (this%has_labels) then
                allocate(this%label(num))
                allocate(this%dilution(num))
                call this%register_int_attribute(this%label, "label", "-")
                call this%register_attribute(this%dilution, "dilution", "-")
            endif

        end subroutine dynamic_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine dynamic_parcel_dealloc(this)
            class(dynamic_parcel_type), intent(inout) :: this

            call try_deallocate(this%volume)
            call try_deallocate(this%vorticity)
            call try_deallocate(this%delta_vor)
            call try_deallocate(this%dilution)
            call try_deallocate(this%label)
            call this%base_dealloc

        end subroutine dynamic_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine dynamic_parcel_resize(this, new_size)
            class(dynamic_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%base_resize(new_size)

            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%vorticity, new_size, this%local_num)
            call resize_array(this%delta_vor, new_size, this%local_num)

            call this%reset_attribute(this%volume, "volume")
            if (this%n_pos == 2) then
                call this%reset_attribute(this%vorticity(1, :), "z_vorticity")
                call this%reset_attribute(this%delta_vor(1, :), "z_vorticity_rk_tendency")
            elseif  (this%n_pos == 3) then
                call this%reset_attribute(this%vorticity(1, :), "x_vorticity")
                call this%reset_attribute(this%vorticity(2, :), "y_vorticity")
                call this%reset_attribute(this%vorticity(3, :), "z_vorticity")
                call this%reset_attribute(this%delta_vor(1, :), "x_vorticity_rk_tendency")
                call this%reset_attribute(this%delta_vor(2, :), "y_vorticity_rk_tendency")
                call this%reset_attribute(this%delta_vor(3, :), "z_vorticity_rk_tendency")
            endif

            if (this%has_labels) then
                call resize_array(this%label, new_size, this%local_num)
                call resize_array(this%dilution, new_size, this%local_num)
                call this%reset_int_attribute(this%label, "label")
                call this%reset_attribute(this%dilution, "dilution")
            endif

        end subroutine dynamic_parcel_resize


        ! Serialize all parcel attributes into a single buffer
        subroutine base_parcel_serialize(this, n, buffer)
            class(base_parcel_type),   intent(in)  :: this
            integer,          intent(in)  :: n
            integer                       :: j
            double precision, intent(out) :: buffer(this%attr_num+this%int_attr_num)

            do j = 1, this%attr_num
                buffer(j)=this%attrib(j)%aptr(n)
            end do

            do j = 1, this%int_attr_num
                buffer(j+this%attr_num)=this%int_attrib(j)%aptr(n)
            end do

        end subroutine base_parcel_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Deserialize parcel attributes into a single buffer
        subroutine base_parcel_deserialize(this, n, buffer)
            class(base_parcel_type),   intent(inout)  :: this
            integer,          intent(in)  :: n
            integer :: j
            double precision, intent(in) :: buffer(this%attr_num+this%int_attr_num)

            do j = 1, this%attr_num
                this%attrib(j)%aptr(n)=buffer(j)
            end do

            do j = 1, this%int_attr_num
                this%int_attrib(j+this%attr_num)%aptr(n)=nint(buffer(j))
            end do

        end subroutine base_parcel_deserialize

   !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_pack(this, pid, num, buffer)
            class(base_parcel_type),   intent(in)  :: this
            integer,          intent(in)  :: pid(:)
            integer,          intent(in)  :: num
            double precision, intent(out) :: buffer(:)
            integer                       :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * (this%attr_num+this%int_attr_num)
                j = n * (this%attr_num+this%int_attr_num)
                call this%serialize(pid(n), buffer(i:j))
            enddo
        end subroutine parcel_pack

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_unpack(this, num, buffer)
            class(base_parcel_type),   intent(inout) :: this
            integer,          intent(in)    :: num
            double precision, intent(in)    :: buffer(:)
            integer                         :: n, i, j

            do n = 1, num
                i = 1 + (n-1) * (this%attr_num+this%int_attr_num)
                j = n * (this%attr_num+this%int_attr_num)
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
            class(base_parcel_type), intent(inout) :: this
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
                print *, "in parcel_container::parcel_delete: more than all parcels are invalid."
                stop
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

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Overwrite parcel n with parcel m. This subroutine only replaces the
        ! common types.
        ! ATTENTION: Extended types must replace additional parcel attributes
        !            in their own routine.
        ! @param[in] n index of parcel to be replaced
        ! @param[in] m index of parcel used to replace parcel at index n
        ! @pre n and m must be valid parcel indices
        subroutine base_parcel_replace(this, n, m)
            class(base_parcel_type),   intent(inout)  :: this
            integer, intent(in) :: n, m
            integer :: j


            do j = 1, this%attr_num
                this%attrib(j)%aptr(n) = this%attrib(j)%aptr(m)
            end do

            do j = 1, this%int_attr_num
                this%int_attrib(j)%aptr(n) = this%int_attrib(j)%aptr(m)
            end do

        end subroutine base_parcel_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Print parcel attributes
        subroutine print_me(this)
            class(base_parcel_type),   intent(inout)  :: this
            integer :: j

            do j = 1, this%attr_num
                print *, this%attrib(j)%name
                print *, this%attrib(j)%unit
                print *, this%attrib(j)%aptr
            end do

            do j = 1, this%int_attr_num
                print *, this%int_attrib(j)%name
                print *, this%int_attrib(j)%unit
                print *, this%int_attrib(j)%aptr
            end do

        end subroutine print_me

end module
