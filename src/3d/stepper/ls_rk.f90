! =============================================================================
!           Low-storage 3rd and 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk
    use options, only : time
    use dimensions, only : I_Z
    use parcel_container
    use parcel_bc
    use parcel_mpi, only : parcel_communicate
    use rk_utils, only: get_dBdt, get_time_step, get_strain_magnitude_field
    use utils, only : write_step
    use parcel_interpl, only : par2grid, grid2par
    use parcel_damping, only : parcel_damp
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use inversion_mod, only : vor2vel, vorticity_tendency
    use parcel_diagnostics, only : calculate_parcel_diagnostics
    use field_diagnostics, only : calculate_field_diagnostics
    use bndry_fluxes, only : apply_bndry_fluxes, bndry_fluxes_time_step
    use mpi_timer, only : start_timer, stop_timer, timings
    use mpi_utils, only : mpi_stop
    implicit none

    private

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: rk_timer

    ! fourth order RK coefficients:
    double precision, dimension(5), target ::             &
        cas4 = (/- 567301805773.0_dp/1357537059087.0_dp,  &
                 -2404267990393.0_dp/2016746695238.0_dp,  &
                 -3550918686646.0_dp/2091501179385.0_dp,  &
                 -1275806237668.0_dp/842570457699.0_dp,   &
                 0.0/) !dummy value, not actually used

    double precision, dimension(5), target ::             &
        cbs4 =  (/1432997174477.0_dp/9575080441755.0_dp,  &
                  5161836677717.0_dp/13612068292357.0_dp, &
                  1720146321549.0_dp/2090206949498.0_dp,  &
                  3134564353537.0_dp/4481467310338.0_dp,  &
                  2277821191437.0_dp/14882151754819.0_dp/)

    ! thrird order RK coefficients:
    double precision, dimension(3), target :: &
        cas3 = (/  -5.0d0 / 9.0d0,            &
                 -153.0d0 / 128.d0,           &
                    0.0 /) !dummy value, not actually used

    double precision, dimension(3), target :: &
        cbs3 =  (/ 1.0d0 / 3.0d0,             &
                  15.0d0 / 16.0d0,            &
                   8.0d0 / 15.0d0 /)

    double precision, dimension(:), pointer :: captr, cbptr

    integer :: n_stages

    public :: ls_rk_setup, rk_timer, ls_rk_step

    contains

        subroutine ls_rk_setup(order)
            integer, intent(in) :: order

            select case (order)
                case (3)
                    captr => cas3
                    cbptr => cbs3
                    n_stages = 3
                case (4)
                    captr => cas4
                    cbptr => cbs4
                    n_stages = 5
                case default
                    call mpi_stop('Only third and fourth order RK supported.')
            end select

        end subroutine ls_rk_setup

        ! Advances the parcels by a single ls-RK-4 step. It calls a
        ! function to obtain the current time step based on the velocity
        ! strain and the buoyancy gradient.
        ! @param[in] t is the time
        ! Precondition: this routine assumes that the fields are
        ! up-to-date for the first sub-step
        subroutine ls_rk_step(t)
            double precision, intent(inout) :: t
            double precision                :: dt
            integer                         :: n

            call par2grid((t > time%initial))

            ! need to be called in order to set initial time step;
            ! this is also needed for the first ls-rk substep
            call vor2vel

            call vorticity_tendency

            ! update the time step
            dt = get_time_step(t)

!            if (dabs(t - time%initial) < 1.0e-13) then
                call bndry_fluxes_time_step(dt)
!            endif

            call grid2par

            call calculate_parcel_diagnostics

            call calculate_field_diagnostics

            call write_step(t)

            call apply_bndry_fluxes(dt)

            do n = 1, n_stages-1
                call ls_rk_substep(dt, n)
            enddo
            call ls_rk_substep(dt, n_stages)

            call start_timer(rk_timer)
            call apply_parcel_reflective_bc
            call stop_timer(rk_timer)

            call get_strain_magnitude_field
            call parcel_damp(dt)

            ! we need to subtract 14 calls since we start and stop
            ! the timer multiple times which increments n_calls
            timings(rk_timer)%n_calls =  timings(rk_timer)%n_calls - (3 * n_stages - 1)

            t = t + dt
        end subroutine ls_rk_step


        ! Do a ls-RK-4 substep.
        ! @param[in] dt is the time step
        ! @param[in] step is the number of the substep (1 to 5 or 1 to 3)
        subroutine ls_rk_substep(dt, step)
            double precision, intent(in) :: dt
            integer,          intent(in) :: step
            double precision             :: ca, cb
            integer                      :: n

            ca = captr(step)
            cb = cbptr(step)

            if (step == 1) then
                call start_timer(rk_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    parcels%delta_b(:, n) = get_dBdt(parcels%B(:, n),           &
                                                     parcels%strain(:, n),      &
                                                     parcels%vorticity(:, n),   &
                                                     parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk_timer)
            else
                call vor2vel

                call vorticity_tendency

                call grid2par(add=.true.)

                call start_timer(rk_timer)

                !$omp parallel do default(shared) private(n)
                do n = 1, n_parcels
                    parcels%delta_b(:, n) = parcels%delta_b(:, n)               &
                                          + get_dBdt(parcels%B(:, n),           &
                                                     parcels%strain(:, n),      &
                                                     parcels%vorticity(:, n),   &
                                                     parcels%volume(n))
                enddo
                !$omp end parallel do

                call stop_timer(rk_timer)
            endif

            call start_timer(rk_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%position(:, n) = parcels%position(:, n) &
                                       + cb * dt * parcels%delta_pos(:, n)

                parcels%vorticity(:, n) = parcels%vorticity(:, n) &
                                        + cb * dt * parcels%delta_vor(:, n)
                parcels%B(:, n) = parcels%B(:, n) &
                                + cb * dt * parcels%delta_b(:, n)
            enddo
            !$omp end parallel do

            call stop_timer(rk_timer)

            if (step == n_stages) then
                call parcel_communicate
               return
            endif

            call start_timer(rk_timer)

            !$omp parallel do default(shared) private(n)
            do n = 1, n_parcels
                parcels%delta_pos(:, n) = ca * parcels%delta_pos(:, n)
                parcels%delta_vor(:, n) = ca * parcels%delta_vor(:, n)
                parcels%delta_b(:, n) = ca * parcels%delta_b(:, n)
            enddo
            !$omp end parallel do

            call parcel_communicate

            call stop_timer(rk_timer)

            call par2grid

        end subroutine ls_rk_substep

end module ls_rk
