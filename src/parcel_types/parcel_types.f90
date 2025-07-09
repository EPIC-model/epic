 module parcel_types
    use physics, only : glat, lambda_c, q_0, qv_dens_coeff, theta_0, gravity
    use constants, only : zero, one
    use parcel_ellipsoid
    implicit none

    type, extends(ellipsoid_parcel_type) :: idealised_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: humidity
        double precision, allocatable, dimension(:) :: buoyancy
        double precision, allocatable, dimension(:) :: merge_humidity
        double precision, allocatable, dimension(:) :: merge_buoyancy

        contains
            procedure :: alloc => idealised_parcel_alloc
            procedure :: dealloc => idealised_parcel_dealloc
            procedure :: resize => idealised_parcel_resize
            procedure :: split => idealised_parcel_split
            procedure :: get_buoyancy => idealised_parcel_get_buoyancy
            procedure :: merge_alloc => idealised_merge_alloc
            procedure :: merge_dealloc => idealised_merge_dealloc
            procedure :: assign_p2m => idealised_assign_p2m
            procedure :: add_p2m => idealised_add_p2m
            procedure :: assign_m2m => idealised_assign_m2m
            procedure :: assign_m2p => idealised_assign_m2p
    end type

    type, extends(ellipsoid_parcel_type) :: realistic_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: qv
        double precision, allocatable, dimension(:) :: ql
        double precision, allocatable, dimension(:) :: theta
        double precision, allocatable, dimension(:) :: Nl ! optional droplet number
        double precision, allocatable, dimension(:) :: merge_qv
        double precision, allocatable, dimension(:) :: merge_ql
        double precision, allocatable, dimension(:) :: merge_theta
        double precision, allocatable, dimension(:) :: merge_Nl

        contains
            procedure :: alloc => realistic_parcel_alloc
            procedure :: dealloc=> realistic_parcel_dealloc
            procedure :: resize => realistic_parcel_resize
            procedure :: split => realistic_parcel_split
            procedure :: get_buoyancy => realistic_parcel_get_buoyancy
            procedure :: merge_alloc => realistic_merge_alloc
            procedure :: merge_dealloc => realistic_merge_dealloc
            procedure :: assign_p2m => realistic_assign_p2m
            procedure :: add_p2m => realistic_add_p2m
            procedure :: assign_m2m => realistic_assign_m2m
            procedure :: assign_m2p => realistic_assign_m2p
    end type

    type, extends(base_parcel_type) :: prec_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: volume
        double precision, allocatable, dimension(:) :: qr
        double precision, allocatable, dimension(:) :: Nr ! droplet number

        contains
            procedure :: alloc => prec_parcel_alloc
            procedure :: dealloc=> prec_parcel_dealloc
            procedure :: resize => prec_parcel_resize

            ! get_buoyancy added here
    end type

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_alloc(this, num)
            class(idealised_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            this%is_idealised = .true.

            call this%ellipsoid_alloc(num)

            allocate(this%buoyancy(num))

            call this%register_attribute(this%buoyancy, "buoyancy", "m/s^2")

            if(this%is_moist) then
                allocate(this%humidity(num))
                call this%register_attribute(this%humidity, "humidity", "kg/kg")
            endif

        end subroutine idealised_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_dealloc(this)
            class(idealised_parcel_type), intent(inout) :: this

            call try_deallocate(this%humidity)
            call try_deallocate(this%buoyancy)
            call this%ellipsoid_dealloc

        end subroutine idealised_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine idealised_parcel_resize(this, new_size)
            class(idealised_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%ellipsoid_resize(new_size)

            call resize_array(this%buoyancy, new_size, this%local_num)
            call this%reset_attribute(this%buoyancy, "buoyancy")

            if(this%is_moist) then
                call resize_array(this%humidity, new_size, this%local_num)
                call this%reset_attribute(this%humidity, "humidity")
            endif

        end subroutine idealised_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_alloc(this, num)
            class(realistic_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            this%is_idealised = .false.

            call this%ellipsoid_alloc(num)

            allocate(this%theta(num))

            call this%register_attribute(this%theta, "theta", "K")

            if(this%is_moist) then
                allocate(this%qv(num))
                allocate(this%ql(num))
                call this%register_attribute(this%qv, "qv", "kg/kg")
                call this%register_attribute(this%ql, "ql", "kg/kg")
            endif

            if(this%has_droplets) then
               allocate(this%Nl(num))
               call this%register_attribute(this%Nl, "Nl", "/m^3")
            endif

        end subroutine realistic_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_dealloc(this)
            class(realistic_parcel_type), intent(inout) :: this

            call try_deallocate(this%theta)

            if(this%is_moist) then
                call try_deallocate(this%qv)
                call try_deallocate(this%ql)
            endif

            if(this%has_droplets) then
                call try_deallocate(this%Nl)
            endif

            call this%ellipsoid_dealloc

        end subroutine realistic_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine realistic_parcel_resize(this, new_size)
            class(realistic_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%ellipsoid_resize(new_size)

            call resize_array(this%theta, new_size, this%local_num)
            call resize_array(this%qv, new_size, this%local_num)
            call resize_array(this%ql, new_size, this%local_num)

            call this%reset_attribute(this%theta, "theta")
            call this%reset_attribute(this%qv, "qv")
            call this%reset_attribute(this%ql, "ql")

            if(this%has_droplets) then
                call resize_array(this%Nl, new_size, this%local_num)
                call this%reset_attribute(this%Nl, "Nl")
            endif

        end subroutine realistic_parcel_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_alloc(this, num)
            class(prec_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: num

            call this%base_alloc(num)

            allocate(this%volume(num))
            allocate(this%qr(num))
            allocate(this%Nr(num))

            call this%register_attribute(this%volume, "volume", "m^3")
            call this%register_attribute(this%qr, "qr", "kg/kg")
            call this%register_attribute(this%Nr, "Nr", "/m^3")

        end subroutine prec_parcel_alloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_dealloc(this)
            class(prec_parcel_type), intent(inout) :: this

            call try_deallocate(this%volume)
            call try_deallocate(this%qr)
            call try_deallocate(this%Nr)

            call this%base_dealloc

        end subroutine prec_parcel_dealloc

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine prec_parcel_resize(this, new_size)
            class(prec_parcel_type), intent(inout) :: this
            integer,        intent(in)    :: new_size

            call this%base_resize(new_size)

            call resize_array(this%volume, new_size, this%local_num)
            call resize_array(this%qr, new_size, this%local_num)
            call resize_array(this%Nr, new_size, this%local_num)

            call this%reset_attribute(this%volume, "volume")
            call this%reset_attribute(this%qr, "qr")
            call this%reset_attribute(this%Nr, "Nr")

        end subroutine prec_parcel_resize

        subroutine realistic_parcel_split(this, n, n_thread_loc, d_pos_split)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            double precision, intent(in) :: d_pos_split(:)

            call this%ellipsoid_split(n, n_thread_loc, d_pos_split)

            this%theta(n_thread_loc) = this%theta(n)
            if(this%is_moist) then
                this%qv(n_thread_loc) = this%qv(n)
                this%ql(n_thread_loc) = this%ql(n)
            endif
            if(this%has_droplets) then
                this%Nl(n_thread_loc) = this%Nl(n)
            endif
        end subroutine realistic_parcel_split

        subroutine idealised_parcel_split(this, n, n_thread_loc,  d_pos_split)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            double precision, intent(in) :: d_pos_split(:)

            call this%ellipsoid_split(n, n_thread_loc, d_pos_split)

            this%buoyancy(n_thread_loc) = this%buoyancy(n)
            if(this%is_moist) then
                this%humidity(n_thread_loc) = this%humidity(n)
            endif
        end subroutine idealised_parcel_split

        pure subroutine realistic_parcel_get_buoyancy(this, num, buoyancy)
            class(realistic_parcel_type), intent(in) :: this
            integer, intent(in) :: num
            double precision, intent(out) :: buoyancy

            if(this%is_moist) then
                ! total buoyancy (including effects of latent heating)
                buoyancy = gravity*(this%theta(num)*(one+qv_dens_coeff*this%qv(num)-this%ql(num))/theta_0-one)
            else
                buoyancy = gravity*(this%theta(num)/theta_0-one)
            endif
        end subroutine realistic_parcel_get_buoyancy

        pure subroutine idealised_parcel_get_buoyancy(this, num, buoyancy)
            class(idealised_parcel_type), intent(in) :: this
            integer, intent(in) :: num
            double precision, intent(out) :: buoyancy
            double precision :: q_c

            !! NOTE: need to check height offset (see previous code)
            if(this%is_moist) then
                q_c = this%humidity(num) &
                    - q_0 * exp(- lambda_c * this%position(this%z_dim, num))
                q_c = max(zero, q_c)
                buoyancy = this%buoyancy(num) + glat * q_c
            else
                buoyancy = this%buoyancy(num)
            endif
        end subroutine idealised_parcel_get_buoyancy

        subroutine realistic_merge_alloc(this, n_merge)
            class(realistic_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: n_merge

            call this%ellipsoid_merge_alloc(n_merge)
            allocate(this%merge_theta(n_merge))
            if(this%is_moist) then
                allocate(this%merge_qv(n_merge))
                allocate(this%merge_ql(n_merge))
            endif
            if(this%has_droplets) then
                allocate(this%merge_Nl(n_merge))
            endif
        end subroutine realistic_merge_alloc

        subroutine idealised_merge_alloc(this, n_merge)
            class(idealised_parcel_type), intent(inout) :: this
            integer,            intent(in)    :: n_merge

            call this%ellipsoid_merge_alloc(n_merge)
            allocate(this%merge_buoyancy(n_merge))
            if(this%is_moist) then
                allocate(this%merge_humidity(n_merge))
            endif

        end subroutine idealised_merge_alloc

        subroutine realistic_merge_dealloc(this)
            class(realistic_parcel_type), intent(inout) :: this

            call this%ellipsoid_merge_dealloc
            call try_deallocate(this%merge_theta)
            if(this%is_moist) then
                call try_deallocate(this%merge_qv)
                call try_deallocate(this%merge_ql)
            endif
            if(this%has_droplets) then
                call try_deallocate(this%merge_Nl)
            endif

        end subroutine realistic_merge_dealloc

        subroutine idealised_merge_dealloc(this)
            class(idealised_parcel_type), intent(inout) :: this

            call this%ellipsoid_merge_dealloc
            call try_deallocate(this%merge_buoyancy)
            if(this%is_moist) then
                call try_deallocate(this%merge_humidity)
            endif

        end subroutine idealised_merge_dealloc

        subroutine idealised_assign_p2m(this, n_in, n_out, temp_volume)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume

            call ellipsoid_assign_p2m(this, n_in, n_out, temp_volume)

            this%merge_buoyancy(n_out) = temp_volume * this%buoyancy(n_in)
            if(this%is_moist) then
                this%merge_humidity(n_out) = temp_volume * this%humidity(n_in)
            endif
        end subroutine idealised_assign_p2m

        subroutine idealised_add_p2m(this, n_in, n_out, temp_volume)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume

            call ellipsoid_add_p2m(this, n_in, n_out, temp_volume)

            this%merge_buoyancy(n_out) = this%merge_buoyancy(n_out) + temp_volume * this%buoyancy(n_in)
            if(this%is_moist) then
                this%merge_humidity(n_out) = this%merge_humidity(n_out) + temp_volume * this%humidity(n_in)
            endif
        end subroutine idealised_add_p2m

        subroutine idealised_assign_m2m(this, n_inout, temp_volume)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_inout
            double precision, intent(in) :: temp_volume

            call ellipsoid_assign_m2m(this, n_inout, temp_volume)

            this%merge_buoyancy(n_inout) = temp_volume * this%merge_buoyancy(n_inout)
            if(this%is_moist) then
                this%merge_humidity(n_inout) = temp_volume * this%merge_humidity(n_inout)
            endif
        end subroutine idealised_assign_m2m

        subroutine idealised_assign_m2p(this, n_in, n_out)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out

            call ellipsoid_assign_m2p(this, n_in, n_out)

            this%buoyancy(n_out) = this%merge_buoyancy(n_in)
            if(this%is_moist) then
                this%humidity(n_out) = this%merge_humidity(n_in)
            endif
        end subroutine idealised_assign_m2p

        subroutine realistic_assign_p2m(this, n_in, n_out, temp_volume)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume

            call ellipsoid_assign_p2m(this, n_in, n_out, temp_volume)

            this%merge_theta(n_out) = temp_volume * this%theta(n_in)
            if(this%is_moist) then
                this%merge_qv(n_out) = temp_volume * this%qv(n_in)
                this%merge_ql(n_out) = temp_volume * this%ql(n_in)
            endif
            if(this%has_droplets) then
                this%merge_Nl(n_out) = temp_volume * this%Nl(n_in)
            endif
        end subroutine realistic_assign_p2m

        subroutine realistic_add_p2m(this, n_in, n_out, temp_volume)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out
            double precision, intent(in) :: temp_volume

            call ellipsoid_add_p2m(this, n_in, n_out, temp_volume)

            this%merge_theta(n_out) = this%merge_theta(n_out) + temp_volume * this%theta(n_in)
            if(this%is_moist) then
                this%merge_qv(n_out) = this%merge_qv(n_out) + temp_volume * this%qv(n_in)
                this%merge_ql(n_out) = this%merge_ql(n_out) + temp_volume * this%ql(n_in)
            endif
            if(this%has_droplets) then
                this%merge_Nl(n_out) = this%merge_Nl(n_out) + temp_volume * this%Nl(n_in)
            endif
        end subroutine realistic_add_p2m

        subroutine realistic_assign_m2m(this, n_inout, temp_volume)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_inout
            double precision, intent(in) :: temp_volume

            call ellipsoid_assign_m2m(this, n_inout, temp_volume)

            this%merge_theta(n_inout) = temp_volume * this%merge_theta(n_inout)
            if(this%is_moist) then
                this%merge_qv(n_inout) = temp_volume * this%merge_qv(n_inout)
                this%merge_ql(n_inout) = temp_volume * this%merge_ql(n_inout)
            endif
            if(this%has_droplets) then
                this%merge_Nl(n_inout) = temp_volume * this%merge_Nl(n_inout)
            endif
        end subroutine realistic_assign_m2m

        subroutine realistic_assign_m2p(this, n_in, n_out)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n_in
            integer, intent(in) :: n_out

            call ellipsoid_assign_m2p(this, n_in, n_out)

            this%theta(n_out) = this%merge_theta(n_in)
            if(this%is_moist) then
                this%qv(n_out) = this%merge_qv(n_in)
                this%ql(n_out) = this%merge_ql(n_in)
            endif
            if(this%has_droplets) then
                this%Nl(n_out) = this%merge_Nl(n_in)
            endif
       end subroutine realistic_assign_m2p


end module
