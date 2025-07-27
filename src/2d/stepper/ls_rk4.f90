! =============================================================================
!               Low-storage 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk4
    use options, only : parcel
    use dynamic_parcels, only : parcels, n_parcels
    use precipitation_parcels, only : prec_parcels, n_prec_parcels
    use parcel_types, only : idealised_parcel_type, realistic_parcel_type
    use parcel_bc
    use rk4_utils, only: get_B, get_time_step
    use utils, only : write_step
    use parcel_interpl, only : par2grid_idealised, par2grid_realistic, grid2par, grid2par_add
    use prec_parcel_interpl, only : prec_par2grid, prec_grid2par, prec_grid2par_add
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg, prec_tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parcel_diagnostics, only : calculate_parcel_diagnostics
    use field_diagnostics, only : calculate_field_diagnostics
    use parameters, only : nx, nz
    use options, only : microphysics
    use timer, only : start_timer, stop_timer, timings
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk4_timer

    double precision, parameter, dimension(5) :: &
        cas = (/- 567301805773.0_dp/1357537059087.0_dp,  &
                -2404267990393.0_dp/2016746695238.0_dp,  &
                -3550918686646.0_dp/2091501179385.0_dp,  &
                -1275806237668.0_dp/842570457699.0_dp,   &
                0.0_dp/) !dummy value, not actually used

    double precision, parameter, dimension(5) :: &
        cbs =  (/1432997174477.0_dp/9575080441755.0_dp,  &
                 5161836677717.0_dp/13612068292357.0_dp, &
                 1720146321549.0_dp/2090206949498.0_dp,  &
                 3134564353537.0_dp/4481467310338.0_dp,  &
                 2277821191437.0_dp/14882151754819.0_dp/)

    contains

        ! Advances the parcels by a single ls-RK-4 step. It calls a
        ! function to obtain the current time step based on the velocity
        ! strain and the buoyancy gradient.
        ! @param[in] t is the time
        ! Precondition: this routine assumes that the fields are
        ! up-to-date for the first sub-step
        subroutine ls_rk4_step(t)
            double precision, intent(inout) :: t
            double precision                :: dt
            integer                         :: n

            select type (parcels)
            type is (idealised_parcel_type)
                call par2grid_idealised(parcels)
            type is (realistic_parcel_type)
                call par2grid_realistic(parcels)
            end select

            if(microphysics%l_precipitation) then
                call prec_par2grid(prec_parcels)
            end if

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            if(microphysics%l_precipitation) then
                call vorticity_tendency(tbuoyg+prec_tbuoyg, vtend)
            else
                call vorticity_tendency(tbuoyg, vtend)
            endif

            ! update the time step
            dt = get_time_step(t)

            call grid2par(parcels%delta_pos, parcels%delta_vor, parcels%strain)

            if(microphysics%l_precipitation) then
                call prec_grid2par(prec_parcels%delta_pos)
                if(microphysics%l_sedimentation) then
                    prec_parcels%local_num = n_prec_parcels
                    call prec_parcels%sedimentation
                endif
            endif

            call calculate_parcel_diagnostics(parcels%delta_pos)

            call calculate_field_diagnostics

            call write_step(t)

            do n = 1, 4
                call ls_rk4_substep(dt, n)

                select type (parcels)
                type is (idealised_parcel_type)
                    call par2grid_idealised(parcels)
                type is (realistic_parcel_type)
                    call par2grid_realistic(parcels)
                end select

                if(microphysics%l_precipitation) then
                    call prec_par2grid(prec_parcels)
                    if(microphysics%l_sedimentation) then
                        prec_parcels%local_num = n_prec_parcels
                        call prec_parcels%sedimentation
                    endif
                end if
            enddo
            call ls_rk4_substep(dt, 5)

            call start_timer(rk4_timer)
            call apply_parcel_bc(parcels%position, parcels%B)
            call stop_timer(rk4_timer)

            ! we need to subtract 14 calls since we start and stop
            ! the timer multiple times which increments n_calls
            timings(rk4_timer)%n_calls =  timings(rk4_timer)%n_calls - 14

            t = t + dt
        end subroutine ls_rk4_step


        ! Do a ls-RK-4 substep.
        ! @param[in] dt is the time step
        ! @param[in] step is the number of the substep (1 to 5)
        subroutine ls_rk4_substep(dt, step)
            double precision, intent(in) :: dt
            integer,          intent(in) :: step
            double precision             :: ca, cb
            integer                      :: n

            ca = cas(step)
            cb = cbs(step)

            if (step == 1) then
                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    parcels%delta_B(:, n) = get_B(parcels%B(:, n), parcels%strain(:, n), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            else
                call vor2vel(vortg, velog, velgradg)

                if(microphysics%l_precipitation) then
                    call vorticity_tendency(tbuoyg+prec_tbuoyg, vtend)
                else
                    call vorticity_tendency(tbuoyg, vtend)
                endif

                call grid2par_add(parcels%delta_pos, parcels%delta_vor, parcels%strain)

                if(microphysics%l_precipitation) then
                    call prec_grid2par_add(prec_parcels%delta_pos)
                endif

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    parcels%delta_b(:, n) = parcels%delta_b(:, n) &
                                 + get_B(parcels%B(:, n), parcels%strain(:, n), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(:, n) = parcels%position(:, n) &
                                      + cb * dt * parcels%delta_pos(:, n)

                parcels%vorticity(1, n) = parcels%vorticity(1, n) + cb * dt * parcels%delta_vor(1, n)
                parcels%B(:, n) = parcels%B(:, n) + cb * dt * parcels%delta_B(:, n)
            enddo
            !$omp end parallel do

            if(microphysics%l_precipitation) then
                !$omp parallel do default(shared) private(n)
                do n = 1, n_prec_parcels
                    prec_parcels%position(:, n) = prec_parcels%position(:, n) &
                                          + cb * dt * prec_parcels%delta_pos(:, n)
                enddo
                !$omp end parallel do
                prec_parcels%local_num = n_prec_parcels
                call prec_parcels%goners
                n_prec_parcels = prec_parcels%local_num
            endif

            call stop_timer(rk4_timer)
            call parcels%saturation_adjustment

            if (step == 5) then
               return
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%delta_pos(:, n) = ca * parcels%delta_pos(:, n)
                parcels%delta_vor(1, n) = ca * parcels%delta_vor(1, n)
                parcels%delta_b(:, n) = ca * parcels%delta_b(:, n)
            enddo
            !$omp end parallel do

            if(microphysics%l_precipitation) then
                !$omp parallel do default(shared) private(n)
                do n = 1, n_prec_parcels
                    prec_parcels%delta_pos(:, n) = ca * prec_parcels%delta_pos(:, n)
                enddo
                !$omp end parallel do
            end if

            call stop_timer(rk4_timer)

        end subroutine ls_rk4_substep

end module ls_rk4
