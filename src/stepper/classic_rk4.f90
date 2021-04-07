! =============================================================================
!               Classical 4th order Runge-Kutta method
!     (see https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
! =============================================================================
module classic_rk4
    use constants, only : max_num_parcels
    use options, only : parcel_info
    use parcel_container
    use parcel_bc
    use ellipse, only : get_B22
    use interpolation, only : grid2par
    use fields, only : strain_f, velocity_f
    implicit none

    type(parcel_container_type) state

    ! classic_rk4 temporaries
    double precision, allocatable, dimension(:, :) :: &
        k1o, k2o, k3o, k4o, &   ! position integration
        strain,             &   ! strain at parcel location
        b1o, b2o, b3o, b4o      ! B matrix integration

    contains

        ! allocate memory of temporaries
        subroutine classic_rk4_alloc(num)
            integer, intent(in) :: num

            allocate(state%position(num, 2))

            allocate(k1o(num, 2))
            allocate(k2o(num, 2))
            allocate(k3o(num, 2))
            allocate(k4o(num, 2))

            if (parcel_info%is_elliptic) then
                allocate(state%B(num, 2))
                allocate(strain(num, 4))

                allocate(b1o(num, 2))
                allocate(b2o(num, 2))
                allocate(b3o(num, 2))
                allocate(b4o(num, 2))
            endif

        end subroutine classic_rk4_alloc

        ! deallocate memory of temporaries
        subroutine classic_rk4_dealloc

            deallocate(state%position)

            deallocate(k1o)
            deallocate(k2o)
            deallocate(k3o)
            deallocate(k4o)

            if (parcel_info%is_elliptic) then
                deallocate(state%B)
                deallocate(strain)

                deallocate(b1o)
                deallocate(b2o)
                deallocate(b3o)
                deallocate(b4o)
            endif

        end subroutine classic_rk4_dealloc


        subroutine classic_rk4_step(dt)
            double precision, intent(in) :: dt

            if (parcel_info%is_elliptic) then
                call classic_rk4_elliptic(dt)
            else
                call classic_rk4_non_elliptic(dt)
            endif

        end subroutine classic_rk4_step


        subroutine classic_rk4_elliptic(dt)
            double precision, intent(in) :: dt

            state%position = parcels%position
            state%B = parcels%B

            call grid2par(state%position, parcels%volume, k1o, velocity_f, state%B)
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B)
            b1o = get_B(state%B, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k1o)


            state%position = parcels%position + 0.5 * dt * k1o
            state%B = parcels%B + 0.5 * dt * b1o

            ! apply position BC
            call apply_parcel_bc(state%position, k1o)

            call grid2par(state%position, parcels%volume, k2o, velocity_f, state%B)
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B)
            b2o = get_B(state%B, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k2o)


            state%position = parcels%position + 0.5 * dt * k2o
            state%B = parcels%B + 0.5 * dt * b2o

            ! apply position BC
            call apply_parcel_bc(state%position, k2o)

            call grid2par(state%position, parcels%volume, k3o, velocity_f, state%B)
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B)
            b3o = get_B(state%B, strain)


            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k3o)

            state%position = parcels%position + dt * k3o
            state%B = parcels%B + dt * b3o

            ! apply position BC
            call apply_parcel_bc(state%position, k3o)

            call grid2par(state%position, parcels%volume, k4o, velocity_f, state%B)
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B)
            b4o = get_B(state%B, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k4o)


            parcels%position = parcels%position &
                             + dt / 6.0 * (k1o + 2.0 * k2o + 2.0 * k3o + k4o)

            parcels%B = parcels%B &
                      + dt / 6.0 * (b1o + 2.0 * b2o + 2.0 * b3o + b4o)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

            ! update parcel velocity
            call grid2par(parcels%position, parcels%volume, &
                          parcels%velocity, velocity_f, parcels%B)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_elliptic


        subroutine classic_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            state%position = parcels%position

            call grid2par(state%position, parcels%volume, k1o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k1o)


            state%position = parcels%position + 0.5 * dt * k1o

            ! apply position BC
            call apply_parcel_bc(state%position, k1o)

            call grid2par(state%position, parcels%volume, k2o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k2o)

            state%position = parcels%position + 0.5 * dt * k2o

            ! apply position BC
            call apply_parcel_bc(state%position, k2o)

            call grid2par(state%position, parcels%volume, k3o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k3o)

            state%position = parcels%position + dt * k3o

            ! apply position BC
            call apply_parcel_bc(state%position, k3o)

            call grid2par(state%position, parcels%volume, k4o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k4o)

            parcels%position = parcels%position &
                             + dt / 6.0 * (k1o + 2.0 * k2o + 2.0 * k3o + k4o)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

            ! update parcel velocity
            call grid2par(parcels%position, parcels%volume, parcels%velocity, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_non_elliptic


        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(max_num_parcels, 2)
            double precision, intent(in) :: S(max_num_parcels, 4)
            double precision             :: Bout(max_num_parcels, 2)
            double precision             :: B22(max_num_parcels)

            B22(1:n_parcels) = get_B22(Bin(1:n_parcels, 1), Bin(1:n_parcels, 2))

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(1:n_parcels, 1) = 2.0 * (S(1:n_parcels, 1) * Bin(1:n_parcels, 1) + &
                                          S(1:n_parcels, 2) * Bin(1:n_parcels, 2))

            ! B12 = dvdx * B11 + dudy * B22
            Bout(1:n_parcels, 2) = S(1:n_parcels, 3) * Bin(1:n_parcels, 1) &
                                 + S(1:n_parcels, 2) * B22(1:n_parcels)

        end function get_B
end module classic_rk4
