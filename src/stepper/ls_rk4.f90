module ls_rk4
    use parameters, only : parcel_info
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B
    use interpolation, only : grid2par, grid2par_add
    use fields, only : strain_f, velocity_f
    implicit none

    double precision, allocatable, dimension(:, :) :: &
        velocity_p, &   ! position integration
        strain, &   ! strain at parcel location
        dbdt      ! B matrix integration

    double precision, parameter :: &
        ca1 = - 567301805773.0/1357537059087.0,  &
        ca2 = -2404267990393.0/2016746695238.0,  &
        ca3 = -3550918686646.0/2091501179385.0,  &
        ca4 = -1275806237668.0/842570457699.0,   &
        ca5 = 0.0,   &  !dummy value, not actually used
        cb1 =  1432997174477.0/9575080441755.0,  &
        cb2 =  5161836677717.0/13612068292357.0, &
        cb3 =  1720146321549.0/2090206949498.0,  &
        cb4 =  3134564353537.0/4481467310338.0,  &
        cb5 =  2277821191437.0/14882151754819.0

    contains

        ! allocate memory of temporaries
        subroutine ls_rk4_alloc(num)
            integer, intent(in) :: num

            allocate(velocity_p(num, 2))

            if (parcel_info%is_elliptic) then
                allocate(strain(num, 4))
                allocate(dbdt(num, 2))
            endif

        end subroutine ls_rk4_alloc

        ! deallocate memory of temporaries
        subroutine ls_rk4_dealloc

            ! TODO

        end subroutine ls_rk4_dealloc


        subroutine ls_rk4_step(dt)
            double precision, intent(in) :: dt

            if (parcel_info%is_elliptic) then
                call ls_rk4_elliptic(dt)
            else
                call ls_rk4_non_elliptic(dt)
            endif

        end subroutine ls_rk4_step


        subroutine ls_rk4_elliptic_substep(ca,cb,dt,l_first,l_last)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            logical, optional, intent(in) :: l_first
            logical, optional, intent(in) :: l_last

            if(present(l_first)) then
                if(l_first) then
                    call grid2par(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
                else
                    call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
                endif
            else
                call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            endif
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            if(present(l_last)) then
               if(l_last) then
                  return
               endif
            endif
            velocity_p(1:n_parcels,:) = ca*velocity_p(1:n_parcels,:)
            dbdt(1:n_parcels,:) = ca*dbdt(1:n_parcels,:)
            return

        end subroutine ls_rk4_elliptic_substep


        subroutine ls_rk4_elliptic(dt)
            double precision, intent(in) :: dt

            call ls_rk4_elliptic_substep(ca1,cb1,dt,l_first=.true.)
            call ls_rk4_elliptic_substep(ca2,cb2,dt)
            call ls_rk4_elliptic_substep(ca3,cb3,dt)
            call ls_rk4_elliptic_substep(ca4,cb4,dt)
            call ls_rk4_elliptic_substep(ca5,cb5,dt,l_last=.true.)

        end subroutine ls_rk4_elliptic


        subroutine ls_rk4_non_elliptic_substep(ca,cb,dt,l_first,l_last)
            double precision, intent(in) :: ca
            double precision, intent(in) :: cb
            double precision, intent(in) :: dt
            logical, optional, intent(in) :: l_first
            logical, optional, intent(in) :: l_last

            if(present(l_first)) then
                if(l_first) then
                    call grid2par(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
                else
                    call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
                endif
            else
                call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            endif
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            if(present(l_last)) then
               if(l_last) then
                  return
               endif
            endif
            velocity_p(1:n_parcels,:) = ca*velocity_p(1:n_parcels,:)
            return

        end subroutine ls_rk4_non_elliptic_substep

        subroutine ls_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            call ls_rk4_non_elliptic_substep(ca1,cb1,dt,l_first=.true.)
            call ls_rk4_non_elliptic_substep(ca2,cb2,dt)
            call ls_rk4_non_elliptic_substep(ca3,cb3,dt)
            call ls_rk4_non_elliptic_substep(ca4,cb4,dt)
            call ls_rk4_non_elliptic_substep(ca5,cb5,dt,l_last=.true.)

        end subroutine ls_rk4_non_elliptic

end module ls_rk4
