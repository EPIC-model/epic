! =============================================================================
!               Low-storage 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk4
    use options, only : parcel
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_dBdt, get_time_step
    use utils, only : write_step
    use parcel_interpl, only : par2grid, grid2par, grid2par_add
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_diagnostics, only : calc_parcel_diagnostics
    use parameters, only : nx, nz
    use timer, only : start_timer, stop_timer, timings
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk4_timer

    double precision, allocatable, dimension(:, :) :: &
        delta_pos, &
        strain,    &   ! strain at parcel location
        delta_b,   &   ! B matrix integration
        delta_vor      ! vorticity integration

    double precision, parameter :: &
        ca1 = - 567301805773.0_dp/1357537059087.0_dp,  &
        ca2 = -2404267990393.0_dp/2016746695238.0_dp,  &
        ca3 = -3550918686646.0_dp/2091501179385.0_dp,  &
        ca4 = -1275806237668.0_dp/842570457699.0_dp,   &
        ca5 = 0.0,   &  !dummy value, not actually used
        cb1 =  1432997174477.0_dp/9575080441755.0_dp,  &
        cb2 =  5161836677717.0_dp/13612068292357.0_dp, &
        cb3 =  1720146321549.0_dp/2090206949498.0_dp,  &
        cb4 =  3134564353537.0_dp/4481467310338.0_dp,  &
        cb5 =  2277821191437.0_dp/14882151754819.0_dp

    contains

        ! Allocate memory of low-storage RK-4 temporaries.
        ! This memory is deallocated at the end of the simulation.
        ! @param[in] num is the size to allocate
        subroutine ls_rk4_alloc(num)
            integer, intent(in) :: num

            allocate(delta_pos(3, num))
            allocate(delta_vor(3, num))
            allocate(strain(5, num))
            allocate(delta_b(5, num))

        end subroutine ls_rk4_alloc

        ! Deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            deallocate(delta_pos)
            deallocate(delta_vor)
            deallocate(strain)
            deallocate(delta_b)

        end subroutine ls_rk4_dealloc

        ! Advances the parcels by a single ls-RK-4 step. It calls a
        ! function to obtain the current time step based on the velocity
        ! strain and the buoyancy gradient.
        ! @param[in] t is the time
        ! Precondition: this routine assumes that the fields are
        ! up-to-date for the first sub-step
        subroutine ls_rk4_step(t)
            double precision, intent(inout) :: t
            double precision                :: dt

            call par2grid

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk4 substep
            call vor2vel(vortg, velog, velgradg)

            call vorticity_tendency(vortg, tbuoyg, velgradg, vtend)

            ! update the time step
            dt = get_time_step(t)

            call grid2par(delta_pos, delta_vor, strain)

            call calc_parcel_diagnostics(delta_pos)

            call write_step(t, dt)

            call ls_rk4_substep(ca1, cb1, dt, 1)

            call par2grid
            call ls_rk4_substep(ca2, cb2, dt, 2)

            call par2grid
            call ls_rk4_substep(ca3, cb3, dt, 3)

            call par2grid
            call ls_rk4_substep(ca4, cb4, dt, 4)

            call par2grid
            call ls_rk4_substep(ca5, cb5, dt, 5)

            call start_timer(rk4_timer)
            call apply_parcel_bc(parcels%position, parcels%B)
            call stop_timer(rk4_timer)

            ! we need to subtract 14 calls since we start and stop
            ! the timer multiple times which increments n_calls
            timings(rk4_timer)%n_calls =  timings(rk4_timer)%n_calls - 14

            t = t + dt
        end subroutine ls_rk4_step


        ! Do a ls-RK-4 substep.
        ! @param[in] ca is a numerical coefficient
        ! @param[in] cb is a numerical coefficient
        ! @param[in] dt is the time step
        ! @param[in] step is the number of the substep (1 to 5)
        subroutine ls_rk4_substep(ca, cb, dt, step)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            integer,          intent(in) :: step
            integer                      :: n


            if (step == 1) then
                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    delta_b(:, n) = get_dBdt(parcels%B(:, n), strain(:, n), &
                                             parcels%vorticity(:, n), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            else
                call vor2vel(vortg, velog, velgradg)

                call vorticity_tendency(vortg, tbuoyg, velgradg, vtend)

                call grid2par_add(delta_pos, delta_vor, strain)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    delta_b(:, n) = delta_b(:, n) &
                                  + get_dBdt(parcels%B(:, n), strain(:, n), &
                                             parcels%vorticity(:, n), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(:, n) = parcels%position(:, n) &
                                       + cb * dt * delta_pos(:, n)

                parcels%vorticity(:, n) = parcels%vorticity(:, n) + cb * dt * delta_vor(:, n)
                parcels%B(:, n) = parcels%B(:, n) + cb * dt * delta_b(:, n)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

            if (step == 5) then
               return
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                delta_pos(:, n) = ca * delta_pos(:, n)
                delta_vor(:, n) = ca * delta_vor(:, n)
                delta_b(:, n) = ca * delta_b(:, n)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

        end subroutine ls_rk4_substep

end module ls_rk4
