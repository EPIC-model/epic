 module parcel_types
    use physics, only : glat, lambda_c, q_0, qv_dens_coeff, theta_0, gravity, r_d, c_p, L_v, p_surf, p_ref, pressure_scale_height
    use constants, only : zero, one
    use timer, only : start_timer, stop_timer
    use parcel_ellipsoid
    use spline_module

    implicit none

    ! For now, put some of the constants for setting up simulations here
    double precision, parameter :: tk0c = 273.15       ! Temperature of freezing in Kelvin
    double precision, parameter :: qsa1 = 3.8          ! Top in equation to calculate qsat
    double precision, parameter :: qsa2 = -17.2693882  ! Constant in qsat equation
    double precision, parameter :: qsa3 = 35.86        ! Constant in qsat equation
    double precision, parameter :: qsa4 = 6.109        ! Constant in qsat equation
    integer :: saturation_adjustment_timer
    logical :: splines_are_initiated = .false.
    type(spline) :: esat_spline, press_spline, exn_spline

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
            procedure :: saturation_adjustment => idealised_saturation_adjustment ! just a stub
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
            procedure :: saturation_adjustment => realistic_saturation_adjustment
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
            procedure :: dealloc => prec_parcel_dealloc
            procedure :: resize => prec_parcel_resize
            procedure :: get_buoyancy => prec_parcel_get_buoyancy

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

       subroutine idealised_saturation_adjustment(this)
            class(idealised_parcel_type), intent(inout) :: this
       end subroutine idealised_saturation_adjustment

       subroutine realistic_saturation_adjustment(this)
            class(realistic_parcel_type), intent(inout) :: this
            double precision :: press, exn, temp, temp_low, qt_start, ql_start, ql_iter
            double precision :: theta_start, temp_start, qsat
            double precision :: err_at_temp, err_at_temp_inv_deriv,efact,divfact, this_height
            double precision :: inv_p_ref, r_d_over_c_p, inv_scale_height, qsat_helper, L_v_over_c_p
            integer :: n, iter, n_parcel_saved

            call start_timer(saturation_adjustment_timer)
            if(.not. this%is_moist) then
                return
            endif

            if(.not. splines_are_initiated) then
                call initiate_thermo_splines
            endif

            inv_scale_height = 1.0d0/pressure_scale_height
            inv_p_ref = 1.0d0/p_ref
            r_d_over_c_p = r_d/c_p
            L_v_over_c_p = L_v/c_p
            n_parcel_saved=this%local_num

            !$omp parallel default(shared)
            !$omp do private(n, press, exn, temp, theta_start, temp_start, ql_start, qt_start, &
            !$omp            temp_low, qsat_helper, efact, qsat, ql_iter, err_at_temp, &
            !$omp            divfact, err_at_temp_inv_deriv, this_height)
            do n = 1, n_parcel_saved
                this_height=this%position(this%n_pos, n)
                press=eval_spline(press_spline, this_height)
                exn=eval_spline(exn_spline, this_height)
                theta_start=this%theta(n)
                temp=theta_start*exn
                temp_start=temp
                ql_start=this%ql(n)
                qt_start=ql_start+this%qv(n)
                ! Test unsaturated case first
                temp_low=temp-L_v_over_c_p*ql_start
                qsat_helper = 0.01d0*press*eval_spline(esat_spline, temp_low) - qsa4
                if(qt_start*qsat_helper < qsa1) then ! Evaporate everything, if needed at all
                   if(ql_start>0.0d0) then
                      this%theta(n)=theta_start-(L_v_over_c_p/exn)*ql_start
                      this%qv(n)=qt_start
                      this%ql(n)=0.0d0
                   end if
                ! Moist case: iterate a few times, start from temp instead of temp_low
                ! Use Newton-Raphson to converge
                else
                   do iter=1,2
                      efact=0.01d0*press*eval_spline(esat_spline, temp)
                      qsat=qsa1/(efact - qsa4)
                      ql_iter=max(qt_start-qsat,0.0d0)
                      err_at_temp=temp-(temp_start-L_v_over_c_p*(ql_start-ql_iter))
                      if(ql_iter>0.0d0) then
                         !calculate 1/(d err/ dt) to save a division latet on
                         divfact=((efact - qsa4)*(efact - qsa4)*(temp - qsa3)*(temp - qsa3))
                         err_at_temp_inv_deriv=divfact/(divfact+L_v_over_c_p*(qsa1*qsa2*efact*(qsa3-tk0c)))
                      else
                         err_at_temp_inv_deriv=1.0d0
                      endif
                      temp=temp-err_at_temp*err_at_temp_inv_deriv
                   enddo
                   qsat=qsa1/(0.01d0*press*eval_spline(esat_spline, temp)  - qsa4)
                   ql_iter=max(qt_start-qsat,0.0d0)
                   this%theta(n)=theta_start-(L_v_over_c_p/exn)*(ql_start-ql_iter)
                   this%qv(n)=qt_start-ql_iter
                   this%ql(n)=ql_iter
                end if
            end do
            !$omp end do
            !$omp end parallel
            call stop_timer(saturation_adjustment_timer)

       end subroutine realistic_saturation_adjustment

       ! can be optimised still
       subroutine set_rh(this, n, rh)
            class(realistic_parcel_type), intent(inout) :: this
            integer, intent(in) :: n
            double precision, intent(in) :: rh
            double precision :: temp, exn, press, qsat

            if(.not. this%is_moist) then
                print *, "Trying to set RH in non-moist simulation"
                stop
            endif

            press=p_surf*exp(-this%position(this%n_pos, n)/pressure_scale_height)
            exn=(press/p_ref)**(r_d/c_p)
            temp=this%theta(n)*exn
            qsat = qsa1/(0.01d0*press*exp(qsa2*(temp - tk0c)/(temp - qsa3)) - qsa4)
            this%qv(n)= rh*qsat

       end subroutine set_rh

  subroutine initiate_thermo_splines
    integer, parameter :: temp_len = 240
    integer, parameter :: height_len = 1000
    double precision :: temperatures(temp_len)
    double precision :: esat_temps(temp_len)
    double precision :: heights(height_len)
    double precision :: press_heights(height_len)
    double precision :: exn_heights(height_len)
    integer :: ii
    double precision :: inv_p_ref, r_d_over_c_p, inv_scale_height, L_v_over_c_p

    inv_scale_height = 1.0d0/pressure_scale_height
    inv_p_ref = 1.0d0/p_ref
    r_d_over_c_p = r_d/c_p
    L_v_over_c_p = L_v/c_p

    do ii=1,temp_len
       temperatures(ii)=160.0d0+ii
       esat_temps(ii)=exp(qsa2*(temperatures(ii) - tk0c)/(temperatures(ii) - qsa3))
    end do
    do ii=1,height_len
       heights(ii)=-40.0d0+ii*40.0d0
       press_heights(ii)=p_surf*exp(-heights(ii)*inv_scale_height)
       exn_heights(ii)=(press_heights(ii)*inv_p_ref)**(r_d_over_c_p)
    end do

    call init_spline(esat_spline, temperatures, esat_temps)
    call init_spline(press_spline, heights, press_heights)
    call init_spline(exn_spline, heights, exn_heights)

  end subroutine initiate_thermo_splines

  pure subroutine prec_parcel_get_buoyancy(this, num, buoyancy)
      class(prec_parcel_type), intent(in) :: this
      integer, intent(in) :: num
      double precision, intent(out) :: buoyancy

      ! precipitation buoyancy term
      buoyancy = -gravity*this%qr(num)

  end subroutine prec_parcel_get_buoyancy

end module
