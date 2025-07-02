 module parcel_ellipsoid

    use parcel_container
    use constants, only : f12
    implicit none

    ! Adding the ellipsoid geomerty to the dynamics
    type, abstract, extends(dynamic_parcel_type) :: ellipsoid_parcel_type ! add procedures for ellipsoid below
        double precision, allocatable, dimension(:,:) :: B
        double precision, allocatable, dimension(:,:) :: delta_B
        double precision, allocatable, dimension(:,:) :: strain
        double precision, allocatable, dimension(:,:) :: Vetas
        double precision, allocatable, dimension(:,:) :: Vtaus
        double precision, allocatable, dimension(:,:) :: evec ! for now add evec for 2D. Need to think about consistency with 3D
        ! Variables for merging. Note these are not attributes as the merges take place on the local domain
        ! Also reallocation is done every time step with length n_merge
        double precision, allocatable, dimension(:,:) :: merge_position
        double precision, allocatable, dimension(:,:) :: merge_B
        double precision, allocatable, dimension(:,:) :: merge_vorticity
        double precision, allocatable, dimension(:) :: merge_volume
        double precision, allocatable, dimension(:) :: merge_dilution
        integer(kind=8), allocatable, dimension(:) :: merge_label

        logical :: par2grid_1p = .false. ! use one point for par2grid
        logical :: grid2par_1p = .false. ! use one point for grid2par
        integer, private :: n_strain = -1      ! number of strain components
        character(len=4), allocatable, dimension(:) :: strain_names ! Names of directions for strain
        integer, private :: n_shape = -1      ! number of shape components
        character(len=3), allocatable, dimension(:) :: shape_names ! Names of shape components
        contains
            procedure :: ellipsoid_alloc => ellipsoid_parcel_alloc
            procedure :: ellipsoid_dealloc => ellipsoid_parcel_dealloc
            procedure :: ellipsoid_resize => ellipsoid_parcel_resize
            procedure :: ellipsoid_split => ellipsoid_parcel_split
            procedure :: ellipsoid_dimensions => set_ellipsoid_dimensions
            procedure :: ellipsoid_merge_alloc => ellipsoid_parcel_merge_alloc
            procedure :: ellipsoid_merge_dealloc => ellipsoid_parcel_merge_dealloc
            procedure(parcel_split), deferred :: split
            procedure(parcel_merge_alloc), deferred :: merge_alloc
            procedure(parcel_merge_dealloc), deferred :: merge_dealloc
            procedure(parcel_assign_p2m), deferred :: assign_p2m
            procedure(parcel_add_p2m), deferred :: add_p2m
            procedure(parcel_assign_m2m), deferred :: assign_m2m
            procedure(parcel_assign_m2p), deferred :: assign_m2p

            ! Other ellipsoid procedures to go here
    end type

    interface
        subroutine parcel_merge_alloc(this, n_merge)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,                 intent(in)    :: n_merge
        end subroutine parcel_merge_alloc

        subroutine parcel_merge_dealloc(this)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
        end subroutine parcel_merge_dealloc

        subroutine parcel_split(this, n, n_thread_loc, d_pos_split)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            double precision, intent(in) :: d_pos_split(:)
        end subroutine parcel_split

        subroutine parcel_assign_p2m(this, n_in, n_out, temp_volume)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume
        end subroutine

        subroutine parcel_add_p2m(this, n_in, n_out, temp_volume)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume
        end subroutine

        subroutine parcel_assign_m2m(this, n_inout, temp_volume)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_inout
            double precision, intent(in) :: temp_volume
        end subroutine

        subroutine parcel_assign_m2p(this,n_in,n_out)
            import ellipsoid_parcel_type
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
        end subroutine

    end interface

    contains

        subroutine set_ellipsoid_dimensions(this)
            class(ellipsoid_parcel_type), intent(inout) :: this

            call this%set_vorticity_dimensions

            if (trim(this%dim_string) == 'xy') then
                this%n_strain=4
                allocate(this%strain_names(this%n_strain))
                this%strain_names(1)='DUDX'
                this%strain_names(2)='DUDY'
                this%strain_names(3)='DVDX'
                this%strain_names(4)='DVDY'
            elseif (trim(this%dim_string) == 'xz') then
                this%n_strain=4
                allocate(this%strain_names(this%n_strain))
                this%strain_names(1)='DUDX'
                this%strain_names(2)='DUDZ'
                this%strain_names(3)='DWDX'
                this%strain_names(4)='DWDZ'
            elseif(trim(this%dim_string) == 'yz') then
                this%n_strain=4
                allocate(this%strain_names(this%n_strain))
                this%strain_names(1)='DVDY'
                this%strain_names(2)='DVDZ'
                this%strain_names(3)='DWDY'
                this%strain_names(4)='DWDZ'
            elseif(trim(this%dim_string) == 'xyz') then
                this%n_strain=9
                allocate(this%strain_names(this%n_strain))
                this%strain_names(1)='DUDX'
                this%strain_names(2)='DUDY'
                this%strain_names(3)='DUDZ'
                this%strain_names(4)='DVDX'
                this%strain_names(5)='DVDY'
                this%strain_names(6)='DVDZ'
                this%strain_names(7)='DWDX'
                this%strain_names(8)='DWDY'
                this%strain_names(9)='DWDZ'
            else
                print *, "ERROR: set_ellipsoid_dimension:: only xy, xz, yz or xyz allowed."
                stop
            endif

            if (len_trim(this%dim_string) == 2) then
                this%n_shape=2
                allocate(this%shape_names(this%n_shape))
                this%shape_names(1)='B11'
                this%shape_names(2)='B12'
            elseif (len_trim(this%dim_string) == 3) then
                this%n_shape=5
                allocate(this%shape_names(this%n_shape))
                this%shape_names(1)='B11'
                this%shape_names(2)='B12'
                this%shape_names(3)='B13'
                this%shape_names(4)='B22'
                this%shape_names(5)='B23'
            else
                print *, "ERROR: set_ellipsoid_dimension:: wrong number of dimensions"
                stop
            endif

        end subroutine set_ellipsoid_dimensions
        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_alloc(this, num)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num
            integer                           :: i_dim

            call this%dynamic_alloc(num)

            allocate(this%B(this%n_shape, num))
            allocate(this%delta_B(this%n_shape, num))
            allocate(this%strain(this%n_strain, num))

            do i_dim=1, this%n_shape
                call this%register_attribute(this%B(i_dim, :), this%shape_names(i_dim), "1/s")
                call this%register_attribute(this%delta_B(i_dim, :), this%shape_names(i_dim)//"_rk_tendency", "1/s")
            enddo

            do i_dim=1, this%n_strain
                call this%register_attribute(this%strain(i_dim, :), this%strain_names(i_dim), "1/s")
            enddo

            if (this%n_pos == 2) then
                allocate(this%evec(2, num))
                call this%register_attribute(this%evec(1, :), "evec 1", "m")
                call this%register_attribute(this%evec(2, :), "evec 2", "m")
            elseif (this%n_pos == 3) then
                allocate(this%Vetas(3, num))
                allocate(this%Vtaus(3, num))
                call this%register_attribute(this%Vetas(1, :), "Veta1", "m")
                call this%register_attribute(this%Vetas(2, :), "Veta2", "m")
                call this%register_attribute(this%Vetas(3, :), "Veta3", "m")
                call this%register_attribute(this%Vtaus(1, :), "Vtau1", "m")
                call this%register_attribute(this%Vtaus(2, :), "Vtau2", "m")
                call this%register_attribute(this%Vtaus(3, :), "Vtau3", "m")
            else
                print *, "ERROR: ellipsoid_parcel_alloc:: n_pos must be 2 or 3."
                stop
            endif

        end subroutine ellipsoid_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_dealloc(this)
            class(ellipsoid_parcel_type), intent(inout) :: this

            call try_deallocate(this%B)
            call try_deallocate(this%delta_B)
            call try_deallocate(this%strain)
            call try_deallocate(this%Vetas)
            call try_deallocate(this%Vtaus)
            call try_deallocate(this%evec)
            call try_deallocate(this%shape_names)
            call try_deallocate(this%strain_names)

            call this%dynamic_dealloc

        end subroutine ellipsoid_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine ellipsoid_parcel_resize(this, new_size)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size
            integer                       :: i_dim

            call this%dynamic_resize(new_size)

            call resize_array(this%B, new_size, this%local_num)
            call resize_array(this%delta_B, new_size, this%local_num)
            call resize_array(this%strain, new_size, this%local_num)

            do i_dim=1,this%n_shape
                call this%reset_attribute(this%B(i_dim, :), this%shape_names(i_dim))
                call this%reset_attribute(this%delta_B(i_dim, :), this%shape_names(i_dim)//"_rk_tendency")
            end do

            do i_dim=1,this%n_strain
                call this%reset_attribute(this%strain(i_dim, :), this%strain_names(i_dim))
            end do

            ! Distinguish between 2D and 3D ellipsoids
            if (this%n_pos == 2) then
                call resize_array(this%evec, new_size, this%local_num)
                call this%reset_attribute(this%evec(1, :), "evec 1")
                call this%reset_attribute(this%evec(2, :), "evec 2")
            elseif (this%n_pos == 3) then
                call resize_array(this%Vetas, new_size, this%local_num)
                call resize_array(this%Vtaus, new_size, this%local_num)
                call this%reset_attribute(this%Vetas(1, :), "Veta1")
                call this%reset_attribute(this%Vetas(2, :), "Veta2")
                call this%reset_attribute(this%Vetas(3, :), "Veta3")
                call this%reset_attribute(this%Vtaus(1, :), "Vtau1")
                call this%reset_attribute(this%Vtaus(2, :), "Vtau2")
                call this%reset_attribute(this%Vtaus(3, :), "Vtau3")
            else
                print *, "ERROR: ellipsoid_parcel_resize:: n_pos must be 2 or 3."
                stop
            endif

        end subroutine ellipsoid_parcel_resize

        subroutine ellipsoid_parcel_split(this, n, n_thread_loc, d_pos_split)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            double precision, intent(in) :: d_pos_split(:)

            this%volume(n) = f12 * this%volume(n)
            this%B(:, n_thread_loc) = this%B(:, n)
            this%vorticity(:, n_thread_loc) = this%vorticity(:, n)
            this%volume(n_thread_loc) = this%volume(n)
            this%position(:, n_thread_loc) = this%position(:, n) - d_pos_split
            this%position(:, n) = this%position(:, n) + d_pos_split
            if(this%has_labels) then
                this%label(n_thread_loc) = this%label(n)
                this%dilution(n_thread_loc) = this%dilution(n)
            endif

        end subroutine ellipsoid_parcel_split

        subroutine ellipsoid_parcel_merge_alloc(this, n_merge)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: n_merge

            allocate(this%merge_position(this%n_pos, n_merge))
            allocate(this%merge_B(this%n_shape+1, n_merge)) ! Extra element for B22/B33
            allocate(this%merge_vorticity(this%n_vor, n_merge))
            allocate(this%merge_volume(n_merge))
            if(this%has_labels) then
                allocate(this%merge_dilution(n_merge))
                allocate(this%merge_label(n_merge))
            endif

        end subroutine ellipsoid_parcel_merge_alloc

        subroutine ellipsoid_parcel_merge_dealloc(this)
            class(ellipsoid_parcel_type), intent(inout) :: this

            call try_deallocate(this%merge_position)
            call try_deallocate(this%merge_B)
            call try_deallocate(this%merge_vorticity)
            call try_deallocate(this%merge_volume)
            if(this%has_labels) then
                call try_deallocate(this%merge_dilution)
                call try_deallocate(this%merge_label)
            endif

        end subroutine ellipsoid_parcel_merge_dealloc

        subroutine ellipsoid_assign_p2m(this, n_in, n_out, temp_volume)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume

            this%merge_vorticity(:, n_out) = temp_volume * this%vorticity(:, n_in)
            if(this%has_labels) then
                this%merge_dilution(n_out) = this%dilution(n_in)
                this%merge_label(n_out) = this%label(n_in)
            endif
        end subroutine ellipsoid_assign_p2m

        subroutine ellipsoid_add_p2m(this, n_in, n_out, temp_volume)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume
            double precision :: rn

            this%merge_vorticity(:, n_out) = this%merge_vorticity(:, n_out) + temp_volume * this%vorticity(:, n_in)
            if(this%has_labels) then
                ! Dilute the parcel when volume is added for now
                ! This could be optimised moving it to later in code
                call random_number(rn)
                if(rn*this%merge_volume(n_out)<temp_volume) then
                   this%merge_dilution(n_out)=this%merge_dilution(n_out)+log(temp_volume/this%merge_volume(n_out))
                else
                   this%merge_label(n_out)=this%label(n_in)
                   this%merge_dilution(n_out)=this%dilution(n_in)+log(temp_volume/this%merge_volume(n_out))
                endif
            endif
        end subroutine ellipsoid_add_p2m

        subroutine ellipsoid_assign_m2m(this, n_inout, temp_volume)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_inout
            double precision, intent(in) :: temp_volume

            this%merge_vorticity(:, n_inout) = temp_volume * this%merge_vorticity(:, n_inout)

        end subroutine ellipsoid_assign_m2m

        subroutine ellipsoid_assign_m2p(this, n_in, n_out)
            class(ellipsoid_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out

            this%vorticity(:, n_out) = this%merge_vorticity(:, n_in)
            if(this%has_labels) then
                this%dilution(n_out) = this%merge_dilution(n_in)
                this%label(n_out) = this%merge_label(n_in)
            endif
       end subroutine ellipsoid_assign_m2p

end module
