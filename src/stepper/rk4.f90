module rk4
    use constants, only : max_num_parcels
    use parameters, only : parcel_info
    use parcel_container
    use parcel_bc
    use ellipse, only : get_B22
    use interpolation, only : grid2par
    use fields, only : strain_f, velocity_f
    implicit none

    type(parcel_container_type) state

    ! rk4 temporaries
    double precision, allocatable, dimension(:, :) :: &
        k1o, k2o, k3o, k4o, &   ! position integration
        s1o, s2o, s3o, s4o, &   ! strain integration
        b1o, b2o, b3o, b4o      ! B matrix integration

    contains

        ! allocate memory of temporaries
        subroutine rk4_alloc(num)
            integer, intent(in) :: num

            allocate(state%position(num, 2))
            allocate(state%B(num, 2))
            allocate(state%volume(num, 1))

            allocate(k1o(num, 2))
            allocate(k2o(num, 2))
            allocate(k3o(num, 2))
            allocate(k4o(num, 2))

            allocate(s1o(num, 4))
            allocate(s2o(num, 4))
            allocate(s3o(num, 4))
            allocate(s4o(num, 4))

            allocate(b1o(num, 2))
            allocate(b2o(num, 2))
            allocate(b3o(num, 2))
            allocate(b4o(num, 2))

        end subroutine rk4_alloc

        ! deallocate memory of temporaries
        subroutine rk4_dealloc

            deallocate(state%position)
            deallocate(state%B)
            deallocate(state%volume)

            deallocate(k1o)
            deallocate(k2o)
            deallocate(k3o)
            deallocate(k4o)

            deallocate(s1o)
            deallocate(s2o)
            deallocate(s3o)
            deallocate(s4o)

            deallocate(b1o)
            deallocate(b2o)
            deallocate(b3o)
            deallocate(b4o)

        end subroutine rk4_dealloc


        subroutine rk4_step(dt)
            double precision, intent(in) :: dt

            if (parcel_info%is_elliptic) then
                call rk4_elliptic(dt)
            else
!                 call rk4_non_elliptic(dt)
            endif

        end subroutine rk4_step


        subroutine rk4_elliptic(dt)
            double precision, intent(in) :: dt

            state%position = parcels%position
            state%B = parcels%B
            state%volume = parcels%volume

            call grid2par(state, k1o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k1o)

            call grid2par(parcels, s1o, strain_f)
            b1o = get_B(parcels%B, s1o)


            state%position = parcels%position + 0.5 * dt + k1o

            ! apply position BC
            call apply_parcel_bc(state%position, k1o)

            call grid2par(state, k2o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k2o)

            call grid2par(parcels, s2o, strain_f)
            b2o = get_B(parcels%B + 0.5 * dt * b1o, s2o)
            state%B = b2o


            state%position = parcels%position + 0.5 * dt * k2o

            ! apply position BC
            call apply_parcel_bc(state%position, k2o)

            print *, "hi 2"
            call grid2par(state, k3o, velocity_f)

            print *, "hi, 2a"

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k3o)

            call grid2par(parcels, s3o, strain_f)
            b3o = get_B(parcels%B + 0.5 * dt * b2o, s3o)
            state%B = b3o

            print *, "hi 3"

            state%position = parcels%position + dt * k3o

            ! apply position BC
            call apply_parcel_bc(state%position, k3o)

            call grid2par(state, k4o, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k4o)

            call grid2par(parcels, s4o, strain_f)
            b4o = get_B(parcels%B + dt * b3o, s4o)

            parcels%position = parcels%position &
                             + dt / 6.0 * (k1o + 2.0 * k2o + 2.0 * k3o + k4o)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

            parcels%B = parcels%B &
                      + dt / 6.0 * (b1o + 2.0 * b2o + 2.0 * b3o + b4o)

            ! update parcel velocity
            call grid2par(parcels, parcels%velocity, velocity_f)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine rk4_elliptic


        subroutine rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

        end subroutine rk4_non_elliptic


        function get_B(Bin, S) result(Bout)
            double precision, intent(in) :: Bin(max_num_parcels, 2)
            double precision, intent(in) :: S(max_num_parcels, 4)
            double precision             :: Bout(max_num_parcels, 2)
            double precision             :: B22(max_num_parcels)

            B22 = get_B22(Bin(:, 1), Bin(:, 2))

            Bout(:, 1) = 2.0 * (S(:, 1) * Bin(:, 1) + S(:, 2) * Bin(:, 2))

            Bout(:, 2) = S(:, 3) * Bin(:, 1) + S(:, 2) * B22

        end function get_B


end module rk4
