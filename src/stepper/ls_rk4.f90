! =============================================================================
!               Low-storage 4th order Runge-Kutta method
!            (see https://doi.org/10.5194/gmd-10-3145-2017)
! =============================================================================
module ls_rk4
    use options, only : parcel
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B, get_stretch
    use parcel_interpl, only : par2grid, grid2par, grid2par_add
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parameters, only : nx, nz
    use timer, only : start_timer, stop_timer
    implicit none

    integer, parameter :: dp=kind(zero)           ! double precision

    integer :: ls_rk4_handle

    double precision, allocatable, dimension(:, :) :: &
        strain, &   ! strain at parcel location
        dbdt        ! B matrix integration

    double precision, allocatable, dimension(:) :: &
        dvordt, &   ! vorticity integration
        dsdt        ! stretch integration (non-elliptic only)

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

            if (parcel%is_elliptic) then
                allocate(dbdt(num, 2))
            else
                allocate(dsdt(num))
            endif

        end subroutine ls_rk4_alloc

        ! deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            ! TODO
            deallocate(dvordt)
            deallocate(strain)

            if (parcel%is_elliptic) then
                deallocate(dbdt)
            else
                deallocate(dsdt)
            endif

        end subroutine ls_rk4_dealloc


        subroutine ls_rk4_step(dt)
            double precision, intent(in) :: dt

            call start_timer(ls_rk4_handle)

            if (parcel%is_elliptic) then
                call ls_rk4_elliptic(dt)
            else
                call ls_rk4_non_elliptic(dt)
            endif

            call stop_timer(ls_rk4_handle)

        end subroutine ls_rk4_step


        subroutine ls_rk4_elliptic_substep(ca, cb, dt, step)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            integer, intent(in) :: step

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)

            if(step==1) then
               call grid2par(parcels%velocity, dvordt, strain)
               dbdt(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                           parcels%volume(1:n_parcels))
            else
               call grid2par_add(parcels%velocity, dvordt, strain)
               dbdt(1:n_parcels,:) = dbdt(1:n_parcels,:) &
                                   + get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                           parcels%volume(1:n_parcels))
            endif
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                                            + cb*dt*parcels%velocity(1:n_parcels,:)
            parcels%vorticity(1:n_parcels) = parcels%vorticity(1:n_parcels) + cb*dt*dvordt(1:n_parcels)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, parcels%velocity)
            if(step==5) then
               return
            endif
            parcels%velocity(1:n_parcels,:) = ca*parcels%velocity(1:n_parcels,:)
            dvordt(1:n_parcels) = ca * dvordt(1:n_parcels)
            dbdt(1:n_parcels,:) = ca*dbdt(1:n_parcels,:)
            return

        end subroutine ls_rk4_elliptic_substep


        subroutine ls_rk4_elliptic(dt)
            double precision, intent(in) :: dt

            call ls_rk4_elliptic_substep(ca1, cb1, dt, 1)
            call ls_rk4_elliptic_substep(ca2, cb2, dt, 2)
            call ls_rk4_elliptic_substep(ca3, cb3, dt, 3)
            call ls_rk4_elliptic_substep(ca4, cb4, dt, 4)
            call ls_rk4_elliptic_substep(ca5, cb5, dt, 5)

        end subroutine ls_rk4_elliptic


        subroutine ls_rk4_non_elliptic_substep(ca, cb, dt, step)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            integer, intent(in) :: step

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)

            if(step==1) then
                call grid2par(parcels%velocity, dvordt, strain)
                dsdt(1:n_parcels) = get_stretch(strain, n_parcels)
            else
                call grid2par_add(parcels%velocity, dvordt, strain)
                dsdt(1:n_parcels) = dsdt(1:n_parcels) + get_stretch(strain, n_parcels)
            endif

            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                                            + cb*dt*parcels%velocity(1:n_parcels,:)
            parcels%vorticity(1:n_parcels) = parcels%vorticity(1:n_parcels) &
                                              + cb*dt*dvordt(1:n_parcels)
            parcels%stretch(1:n_parcels) = parcels%stretch(1:n_parcels) &
                                         + cb * dt * dsdt(1:n_parcels)

            call apply_parcel_bc(parcels%position, parcels%velocity)

            if(step==5) then
               return
            endif
            parcels%velocity(1:n_parcels,:) = ca*parcels%velocity(1:n_parcels,:)
            dvordt(1:n_parcels) = ca * dvordt(1:n_parcels)
            dsdt(1:n_parcels) = ca * dsdt(1:n_parcels)
            return

        end subroutine ls_rk4_non_elliptic_substep

        subroutine ls_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            call ls_rk4_non_elliptic_substep(ca1, cb1, dt, 1)
            call ls_rk4_non_elliptic_substep(ca2, cb2, dt, 2)
            call ls_rk4_non_elliptic_substep(ca3, cb3, dt, 4)
            call ls_rk4_non_elliptic_substep(ca4, cb4, dt, 4)
            call ls_rk4_non_elliptic_substep(ca5, cb5, dt, 5)

        end subroutine ls_rk4_non_elliptic

end module ls_rk4
