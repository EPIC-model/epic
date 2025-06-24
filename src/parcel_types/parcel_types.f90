 module parcel_types

    use parcel_ellipsoid
    implicit none

    type, extends(ellipsoid_parcel_type) :: idealised_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: humidity
        double precision, allocatable, dimension(:) :: buoyancy
        logical :: is_moist = .false.

        contains
            procedure :: alloc => idealised_parcel_alloc
            procedure :: dealloc=> idealised_parcel_dealloc
            procedure :: resize => idealised_parcel_resize
            procedure :: split => idealised_parcel_split
            procedure :: get_buoyancy => idealised_parcel_get_buoyancy

    end type

    type, extends(ellipsoid_parcel_type) :: realistic_parcel_type ! add procedures
        double precision, allocatable, dimension(:) :: qv
        double precision, allocatable, dimension(:) :: ql
        double precision, allocatable, dimension(:) :: theta
        double precision, allocatable, dimension(:) :: Nl ! optional droplet number
        logical :: is_moist = .false.
        logical :: has_droplets = .false.

        contains
            procedure :: alloc => realistic_parcel_alloc
            procedure :: dealloc=> realistic_parcel_dealloc
            procedure :: resize => realistic_parcel_resize
            procedure :: split => realistic_parcel_split
            procedure :: get_buoyancy => realistic_parcel_get_buoyancy

            ! get_buoyancy added here
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
            call try_deallocate(this%qv)
            call try_deallocate(this%ql)
            call try_deallocate(this%Nl)

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

        subroutine realistic_parcel_split(this, n, n_thread_loc)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            call this%ellipsoid_split(n, n_thread_loc)
        end subroutine realistic_parcel_split

        subroutine idealised_parcel_split(this, n, n_thread_loc)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            integer, intent(in) :: n_thread_loc
            call this%ellipsoid_split(n, n_thread_loc)
        end subroutine idealised_parcel_split

        pure subroutine realistic_parcel_get_buoyancy(this, num, buoyancy)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: num
            double precision, intent(out) :: buoyancy
            buoyancy = 0.0
        end subroutine realistic_parcel_get_buoyancy

        pure subroutine idealised_parcel_get_buoyancy(this, num, buoyancy)
            class(idealised_parcel_type), intent(inout) :: this
            integer, intent(in) :: num
            double precision, intent(out) :: buoyancy
            buoyancy = 0.0
        end subroutine idealised_parcel_get_buoyancy

end module
