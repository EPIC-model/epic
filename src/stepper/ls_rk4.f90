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
    use parameters, only : nx, nz
    use timer, only : start_timer, stop_timer, timings
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk4_timer

    double precision, allocatable, dimension(:, :) :: &
        strain, &   ! strain at parcel location
        dbdt        ! B matrix integration

    double precision, allocatable, dimension(:) :: &
        dvordt      ! vorticity integration

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

            allocate(dvordt(num))
            allocate(strain(num, 4))
            allocate(dbdt(num, 2))

        end subroutine ls_rk4_alloc

        ! deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            deallocate(dvordt)
            deallocate(strain)
            deallocate(dbdt)

        end subroutine ls_rk4_dealloc


        subroutine ls_rk4_step(dt)
            double precision, intent(in) :: dt

            call ls_rk4_substep(ca1, cb1, dt, 1)
            call ls_rk4_substep(ca2, cb2, dt, 2)
            call ls_rk4_substep(ca3, cb3, dt, 3)
            call ls_rk4_substep(ca4, cb4, dt, 4)
            call ls_rk4_substep(ca5, cb5, dt, 5)

            call start_timer(rk4_timer)

            call apply_all_reflective_bc(parcels%position(1:n_parcels, :), &
                                         parcels%B(1:n_parcels, :))

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

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)

            if (step == 1) then
                call grid2par(parcels%velocity, dvordt, strain)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    dbdt(n,:) = get_B(parcels%B(n,:), strain(n,:), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            else
                call grid2par_add(parcels%velocity, dvordt, strain)

                call start_timer(rk4_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    dbdt(n,:) = dbdt(n,:) &
                              + get_B(parcels%B(n,:), strain(n,:), parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk4_timer)
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(n,:) = parcels%position(n,:) &
                                      + cb * dt * parcels%velocity(n,:)

                call apply_periodic_bc(parcels%position(n,:))

                parcels%vorticity(n) = parcels%vorticity(n) + cb * dt * dvordt(n)
                parcels%B(n,:) = parcels%B(n,:) + cb * dt * dbdt(n,:)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

            if (step == 5) then
               return
            endif

            call start_timer(rk4_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%velocity(n,:) = ca * parcels%velocity(n,:)
                dvordt(n) = ca * dvordt(n)
                dbdt(n,:) = ca * dbdt(n,:)
            enddo
            !$omp end parallel do

            call stop_timer(rk4_timer)

        end subroutine ls_rk4_substep

end module ls_rk4
