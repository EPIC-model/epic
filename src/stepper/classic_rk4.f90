! =============================================================================
!               Classical 4th order Runge-Kutta method
!     (see https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
! =============================================================================
module classic_rk4
    use constants, only : max_num_parcels
    use options, only : parcel_info
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B
    use interpolation, only : grid2par
    use fields, only : strain_f, velocity_f
    implicit none
    integer, parameter :: dp=kind(0.d0)           ! double precision

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

            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:)
            state%B(1:n_parcels,:) = parcels%B(1:n_parcels,:)

            call grid2par(state%position, parcels%volume, k1o, velocity_f, state%B, exact='velocity')
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B, exact='strain')
            b1o(1:n_parcels,:) = get_B(state%B(1:n_parcels,:), strain(1:n_parcels,:))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k1o)


            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + 0.5_dp * dt * k1o(1:n_parcels,:)
            state%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + 0.5_dp * dt * b1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k1o)

            call grid2par(state%position, parcels%volume, k2o, velocity_f, state%B, exact='velocity')
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B, exact='strain')
            b2o(1:n_parcels,:) = get_B(state%B(1:n_parcels,:), strain(1:n_parcels,:))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k2o)


            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + 0.5_dp * dt * k2o(1:n_parcels,:)
            state%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + 0.5_dp * dt * b2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k2o)

            call grid2par(state%position, parcels%volume, k3o, velocity_f, state%B, exact='velocity')
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B, exact='strain')
            b3o(1:n_parcels,:) = get_B(state%B(1:n_parcels,:), strain(1:n_parcels,:))


            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k3o)

            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) + dt * k3o(1:n_parcels,:)
            state%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) + dt * b3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k3o)

            call grid2par(state%position, parcels%volume, k4o, velocity_f, state%B, exact='velocity')
            call grid2par(state%position, parcels%volume, strain, strain_f, state%B, exact='strain')
            b4o(1:n_parcels,:) = get_B(state%B(1:n_parcels,:), strain(1:n_parcels,:))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k4o)


            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                             + dt / 6.0_dp * (k1o(1:n_parcels,:) + 2.0_dp &
                             * k2o(1:n_parcels,:) + 2.0_dp * k3o(1:n_parcels,:) + k4o(1:n_parcels,:))

            parcels%B(1:n_parcels,:) = parcels%B(1:n_parcels,:) &
                      + dt / 6.0_dp * (b1o(1:n_parcels,:) + 2.0_dp * b2o(1:n_parcels,:) &
                      + 2.0_dp * b3o(1:n_parcels,:) + b4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

            ! update parcel velocity
            call grid2par(parcels%position, parcels%volume, &
                          parcels%velocity, velocity_f, parcels%B, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_elliptic


        subroutine classic_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            state%position = parcels%position

            call grid2par(state%position, parcels%volume, k1o, velocity_f, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k1o)


            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                                          + 0.5_dp * dt * k1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k1o)

            call grid2par(state%position, parcels%volume, k2o, velocity_f, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k2o)

            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                                          + 0.5_dp * dt * k2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k2o)

            call grid2par(state%position, parcels%volume, k3o, velocity_f, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k3o)

            state%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                                            + dt * k3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(state%position, k3o)

            call grid2par(state%position, parcels%volume, k4o, velocity_f, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(state%position, k4o)

            parcels%position(1:n_parcels,:) = parcels%position(1:n_parcels,:) &
                             + dt / 6.0_dp * (k1o(1:n_parcels,:) + 2.0_dp &
                             * k2o(1:n_parcels,:) + 2.0_dp * k3o(1:n_parcels,:) + k4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

            ! update parcel velocity
            call grid2par(parcels%position, parcels%volume, parcels%velocity, velocity_f, exact='velocity')

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_non_elliptic

end module classic_rk4
