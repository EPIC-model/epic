module ls_rk4
    use constants, only : max_num_parcels
    use parameters, only : parcel_info
    use parcel_container
    use parcel_bc
    use ellipse, only : get_B22
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


        subroutine ls_rk4_elliptic(dt)
            double precision, intent(in) :: dt

            call grid2par(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb1*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb1*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca1*velocity_p(1:n_parcels,:)
            dbdt(1:n_parcels,:) = ca1*dbdt(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = dbdt+get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb2*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb2*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca2*velocity_p(1:n_parcels,:)
            dbdt(1:n_parcels,:) = ca2*dbdt(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = dbdt+get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb3*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb3*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca3*velocity_p(1:n_parcels,:)
            dbdt(1:n_parcels,:) = ca3*dbdt(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = dbdt+get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb4*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb4*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca4*velocity_p(1:n_parcels,:)
            dbdt(1:n_parcels,:) = ca4*dbdt(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f, parcels%B)
            call grid2par(parcels%position, parcels%volume, strain, strain_f, parcels%B)
            dbdt = dbdt+get_B(parcels%B, strain)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb5*dt*velocity_p(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + cb5*dt*dbdt(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)

        end subroutine ls_rk4_elliptic


        subroutine ls_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            call grid2par(parcels%position, parcels%volume, velocity_p, velocity_f)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb1*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca1*velocity_p(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb2*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca2*velocity_p(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb3*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca3*velocity_p(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb4*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)
            velocity_p(1:n_parcels,:) = ca4*velocity_p(1:n_parcels,:)

            call grid2par_add(parcels%position, parcels%volume, velocity_p, velocity_f)
            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + cb5*dt*velocity_p(1:n_parcels,:)
            call apply_parcel_bc(parcels%position, velocity_p)

        end subroutine ls_rk4_non_elliptic

        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(max_num_parcels, 2)
            double precision, intent(in) :: S(max_num_parcels, 4)
            double precision             :: Bout(max_num_parcels, 2)
            double precision             :: B22(max_num_parcels)

            B22 = get_B22(Bin(:, 1), Bin(:, 2))

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(:, 1) = 2.0 * (S(:, 1) * Bin(:, 1) + S(:, 2) * Bin(:, 2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(:, 2) = S(:, 3) * Bin(:, 1) + S(:, 2) * B22

        end function get_B

end module ls_rk4
