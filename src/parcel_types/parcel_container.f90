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
        integer             :: n_pos = -1   ! number of spatial dimensions
        character(len=1), allocatable, dimension(:) :: pos_names ! Names of directions for positions
        character(len=8) :: dim_string ! Names of directions for positions
        integer :: x_dim = -1 ! For reverse lookup
        integer :: y_dim = -1 ! For reverse lookup
        integer :: z_dim = -1 ! For reverse lookup

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
            procedure :: set_dimensions
            procedure :: register_attribute
            procedure :: register_int_attribute
            procedure :: reset_attribute
            procedure :: reset_int_attribute
    end type

    ! This type is for parcels associated with dynamics (which will always have volume and voriticity in EPIC)
    type, abstract, extends(base_parcel_type) :: dynamic_parcel_type ! add procedures
        double precision, allocatable, dimension(:)   :: volume
        double precision, allocatable, dimension(:,:) :: vorticity
        double precision, allocatable, dimension(:,:) :: delta_vor
        double precision, allocatable, dimension(:) :: dilution
        integer(kind=8), allocatable, dimension(:) :: label
        integer :: n_vor = -1      ! number of voriticity components
        character(len=1), allocatable, dimension(:) :: vor_names ! Names of vorticity components
        logical   :: has_labels = .false.
        logical   :: is_idealised = .false.
        logical :: is_moist = .false.
        logical :: has_droplets = .false.

        contains
            procedure :: dynamic_alloc => dynamic_parcel_alloc
            procedure :: dynamic_dealloc => dynamic_parcel_dealloc
            procedure :: dynamic_resize => dynamic_parcel_resize
            procedure :: set_vorticity_dimensions
            procedure(dynamic_parcel_get_buoyancy), deferred :: get_buoyancy
            procedure(dynamic_parcel_saturation_adjustment), deferred :: saturation_adjustment

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
        module procedure :: try_deallocate_string
    end interface try_deallocate

    interface
        subroutine dynamic_parcel_get_buoyancy(this, num, buoyancy)
            import dynamic_parcel_type
            class(dynamic_parcel_type), intent(in) :: this
            integer,                 intent(in)    :: num
            double precision, intent(out) :: buoyancy
        end subroutine dynamic_parcel_get_buoyancy
    end interface

    interface
        subroutine dynamic_parcel_saturation_adjustment(this)
            import dynamic_parcel_type
            class(dynamic_parcel_type), intent(inout) :: this
        end subroutine dynamic_parcel_saturation_adjustment
    end interface

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

        subroutine try_deallocate_string(in_array)
            character(len=*), allocatable, dimension(:) :: in_array

            if (allocated(in_array)) then
                deallocate(in_array)
            endif

        end subroutine try_deallocate_string

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Allocate parcel memory
        ! ATTENTION: Extended types must allocate additional parcel attributes
        !            in their own routine.
        ! @param[in] num number of parcels
        subroutine base_parcel_alloc(this, num)
            class(base_parcel_type), intent(inout), target :: this
            integer,        intent(in)    :: num
            integer                       :: n

            this%max_num = num
            this%local_num = num ! Usually this is done in parcel_init
            this%attr_num = 0
            this%int_attr_num = 0

            if (this%n_pos > 3) then
                print *, "Only 3 dimensions allowed."
                stop
            endif

            allocate(this%position(this%n_pos, num))
            allocate(this%delta_pos(this%n_pos, num))

            do n = 1, this%n_pos
                call this%register_attribute(this%position(n, :), this%pos_names(n) // "_position", "m")
                call this%register_attribute(this%delta_pos(n, :), this%pos_names(n) // "_position_rk_tendency", "m/s")
            enddo

            do n = 1, this%n_pos
               if(this%pos_names(n)=='x') then
                  this%x_dim=n
               elseif(this%pos_names(n)=='y') then
                  this%y_dim=n
               elseif(this%pos_names(n)=='z') then
                  this%z_dim=n
               endif
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
            call try_deallocate(this%pos_names)

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
            integer                       :: n

            if (new_size < this%local_num) then
                print *, "in parcel_container::base_parcel_resize: losing parcels when resizing."
                stop
            endif

            this%max_num = new_size

            if ((this%n_pos < 1) .or. (this%n_pos > 3)) then
                print *, "ERROR: base_parcel_resize:: only 1-, 2- or 3-dimensional parcels allowed."
                stop
            endif

            call resize_array(this%position, new_size, this%local_num)
            call resize_array(this%delta_pos, new_size, this%local_num)

            do n = 1, this%n_pos
                call this%reset_attribute(this%position(n, :), this%pos_names(n) // "_position")
                call this%reset_attribute(this%delta_pos(n, :), this%pos_names(n) // "_position_rk_tendency")
            enddo

        end subroutine base_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_dimensions(this)
            class(base_parcel_type), intent(inout) :: this
            integer :: n_dim, i_dim

            n_dim = len_trim(this%dim_string)

            if (this%n_pos /= -1) then
                print *, "WARNING: Dimension already set."
                return
            endif

            if ((n_dim < 1) .or. (n_dim > 3)) then
                print *, "ERROR: set_dimension:: only 1-, 2- or 3-dimensional parcels allowed."
                stop
            endif

            this%n_pos=n_dim

            allocate(this%pos_names(this%n_pos))

            do i_dim = 1,this%n_pos
               this%pos_names(i_dim)=this%dim_string(i_dim:i_dim)
            enddo

        end subroutine set_dimensions

        subroutine set_vorticity_dimensions(this)
            class(dynamic_parcel_type), intent(inout) :: this
            integer :: i_dim

            call set_dimensions(this)

            if (trim(this%dim_string) == trim('xy')) then
                this%n_vor=1
                allocate(this%vor_names(this%n_vor))
                this%vor_names(1)='z'
            elseif (trim(this%dim_string) == trim('xz')) then
                this%n_vor=1
                allocate(this%vor_names(this%n_vor))
                this%vor_names(1)='y'
            elseif(trim(this%dim_string) == trim('yz')) then
                this%n_vor=1
                allocate(this%vor_names(this%n_vor))
                this%vor_names(1)='x'
            elseif(trim(this%dim_string) == trim('xyz')) then
                this%n_vor=3
                allocate(this%vor_names(this%n_vor))
                do i_dim = 1,this%n_vor
                   this%vor_names(i_dim)=this%dim_string(i_dim:i_dim)
                enddo
            else
                print *, "ERROR: set_vorticiy_dimension:: only xy, xz, yz or xyz allowed."
                stop
            endif

        end subroutine set_vorticity_dimensions

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
                     print *, name
                     stop
                 end if
            end do

            do j = 1, this%int_attr_num
                 if(trim(this%int_attrib(j)%name) == trim(name)) then
                     print *, "Attribute name not unique"
                     print *, name
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
            integer                           :: i_dim

            call this%base_alloc(num)

            allocate(this%volume(num))
            call this%register_attribute(this%volume, "volume", "m^3")

            if(.not. ((this%n_vor == 1) .or. (this%n_vor == 3))) then
                print *, "Vorticity needs 1 or 3 components."
                stop
            endif

            allocate(this%vorticity(this%n_vor, num))
            allocate(this%delta_vor(this%n_vor, num))
            do i_dim=1,this%n_vor
                call this%register_attribute(this%vorticity(i_dim, :), this%vor_names(i_dim)//"_vorticity", "1/s")
                call this%register_attribute(this%delta_vor(i_dim, :), this%vor_names(i_dim)//"_vorticity_rk_tendency", "1/s")
            enddo

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
            if(this%has_labels) then
                call try_deallocate(this%dilution)
                call try_deallocate(this%label)
            endif
            call try_deallocate(this%vor_names)
            call this%base_dealloc

        end subroutine dynamic_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine dynamic_parcel_resize(this, new_size)
            class(dynamic_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size
            integer                       :: i_dim


            call this%base_resize(new_size)

            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%vorticity, new_size, this%local_num)
            call resize_array(this%delta_vor, new_size, this%local_num)

            call this%reset_attribute(this%volume, "volume")

            do i_dim=1,this%n_vor
                call this%reset_attribute(this%vorticity(i_dim, :), this%vor_names(i_dim)//"_vorticity")
                call this%reset_attribute(this%delta_vor(i_dim, :), this%vor_names(i_dim)//"_vorticity_rk_tendency")
            end do

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
