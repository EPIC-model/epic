! =============================================================================
!               Low-storage 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk4
    use options, only : parcel
    use parcel_container
    use surface_parcel_container, only : n_top_parcels, top_parcels &
                                       , n_bot_parcels, bot_parcels
    use parcel_bc
    use rk4_utils, only: get_B, get_time_step
    use utils, only : write_step
    use parcel_interpl, only : par2grid, grid2par
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parcel_diagnostics, only : calculate_parcel_diagnostics
    use field_diagnostics, only : calculate_field_diagnostics
    use parameters, only : nx, nz
    use timer, only : start_timer, stop_timer, timings
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk4_timer

    double precision, parameter, dimension(5) :: &
        cas = (/- 567301805773.0_dp/1357537059087.0_dp,  &
                -2404267990393.0_dp/2016746695238.0_dp,  &
                -3550918686646.0_dp/2091501179385.0_dp,  &
                -1275806237668.0_dp/842570457699.0_dp,   &
                0.0/) !dummy value, not actually used

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

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(tbuoyg, vtend)

            ! update the time step
            dt = get_time_step(t)

            call grid2par

            call calculate_parcel_diagnostics(parcels%delta_pos)

            call calculate_field_diagnostics

            call write_step(t)

            do n = 1, 4
                call ls_rk4_substep(dt, n)
                call apply_surface_parcel_bc
                call par2grid
!                 call write_step(t+dble(n)*0.001d0)
            enddo
            call ls_rk4_substep(dt, 5)
!             call write_step(t+0.005d0)

            call start_timer(rk4_timer)
            call apply_parcel_bc
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
                    parcels%delta_b(:, n) = get_B(parcels%B(:, n),      &
                                                  parcels%strain(:, n), &
                                                  parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            else
                call vor2vel(vortg, velog, velgradg)

                call vorticity_tendency(tbuoyg, vtend)

                call grid2par(add=.true.)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    parcels%delta_b(:, n) = parcels%delta_b(:, n)       &
                                          + get_B(parcels%B(:, n),      &
                                                  parcels%strain(:, n), &
                                                  parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(:, n) = parcels%position(:, n) &
                                       + cb * dt * parcels%delta_pos(:, n)

                parcels%vorticity(n) = parcels%vorticity(n) &
                                     + cb * dt * parcels%delta_vor(n)
                parcels%B(:, n) = parcels%B(:, n) &
                                + cb * dt * parcels%delta_b(:, n)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(n)
            do n = 1, n_bot_parcels
                bot_parcels%position(n) = bot_parcels%position(n) &
                                        + cb * dt * bot_parcels%delta_pos(n)
                bot_parcels%vorticity(n) = bot_parcels%vorticity(n) &
                                         + cb * dt * bot_parcels%delta_vor(n)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(n)
            do n = 1, n_top_parcels
                top_parcels%position(n) = top_parcels%position(n) &
                                        + cb * dt * top_parcels%delta_pos(n)
                top_parcels%vorticity(n) = top_parcels%vorticity(n) &
                                         + cb * dt * top_parcels%delta_vor(n)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

            if (step == 5) then
               return
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%delta_pos(:, n) = ca * parcels%delta_pos(:, n)
                parcels%delta_vor(n) = ca * parcels%delta_vor(n)
                parcels%delta_b(:, n) = ca * parcels%delta_b(:, n)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(n)
            do n = 1, n_bot_parcels
                bot_parcels%delta_pos(n) = ca * bot_parcels%delta_pos(n)
                bot_parcels%delta_vor(n) = ca * bot_parcels%delta_vor(n)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(n)
            do n = 1, n_top_parcels
                top_parcels%delta_pos(n) = ca * top_parcels%delta_pos(n)
                top_parcels%delta_vor(n) = ca * top_parcels%delta_vor(n)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

        end subroutine ls_rk4_substep

end module ls_rk4
