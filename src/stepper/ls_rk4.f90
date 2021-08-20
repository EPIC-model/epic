! =============================================================================
!               Low-storage 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk4
    use options, only : parcel
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B
    use parcel_interpl, only : par2grid, grid2par, grid2par_add
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parcel_diagnostics, only : calc_parcel_diagnostics
    use parameters, only : nx, nz
    use timer, only : start_timer, stop_timer, timings
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk4_timer

    double precision, allocatable, dimension(:, :) :: &
        delta_pos, &
        strain,    &   ! strain at parcel location
        delta_b        ! B matrix integration

    double precision, allocatable, dimension(:) :: &
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

        ! allocate memory of temporaries
        subroutine ls_rk4_alloc(num)
            integer, intent(in) :: num

            allocate(delta_pos(num, 2))
            allocate(delta_vor(num))
            allocate(strain(num, 4))
            allocate(delta_b(num, 2))

        end subroutine ls_rk4_alloc

        ! deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            deallocate(delta_pos)
            deallocate(delta_vor)
            deallocate(strain)
            deallocate(delta_b)

        end subroutine ls_rk4_dealloc

        ! Precondition: this routine assumes that the fields are
        ! up-to-date for the first sub-step
        subroutine ls_rk4_step(dt)
            double precision, intent(in) :: dt

            ! no need to call par2grid (fields are up-to-date)
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

        end subroutine ls_rk4_step


        subroutine ls_rk4_substep(ca, cb, dt, step)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            integer,          intent(in) :: step
            integer                      :: n

            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)

            if (step == 1) then
                call grid2par(delta_pos, delta_vor, strain)

                call calc_parcel_diagnostics(delta_pos)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    delta_b(n,:) = get_B(parcels%B(n,:), strain(n,:), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            else
                call grid2par_add(delta_pos, delta_vor, strain)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    delta_b(n,:) = delta_b(n,:) &
                                 + get_B(parcels%B(n,:), strain(n,:), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(n,:) = parcels%position(n,:) &
                                      + cb * dt * delta_pos(n,:)

                parcels%vorticity(n) = parcels%vorticity(n) + cb * dt * delta_vor(n)
                parcels%B(n,:) = parcels%B(n,:) + cb * dt * delta_b(n,:)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

            if (step == 5) then
               return
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                delta_pos(n,:) = ca * delta_pos(n,:)
                delta_vor(n) = ca * delta_vor(n)
                delta_b(n,:) = ca * delta_b(n,:)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

        end subroutine ls_rk4_substep

end module ls_rk4
