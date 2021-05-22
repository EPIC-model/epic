! =============================================================================
!               Classical 4th order Runge-Kutta method
!     (see https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
! =============================================================================
module classic_rk4
    use constants, only : max_num_parcels
    use parameters, only : nx, nz
    use options, only : parcel_info
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B
    use tri_inversion, only : vor2vel
    use parcel_interpl, only : par2grid, grid2par
    use fields, only : velgradg, velog, vortg, vtend
    implicit none
    integer, parameter :: dp=kind(0.d0)           ! double precision

    ! classic_rk4 temporaries
    double precision, allocatable, dimension(:, :) :: &
        k1o, k2o, k3o, k4o,     &   ! position integration
        w1o, w2o, w3o, w4o,     &   ! vorticity integration
        strain,                 &   ! strain at parcel location
        b1o, b2o, b3o, b4o,     &   ! B matrix integration
        inipos, inivor, iniB        ! input parcel attributes (before RK4 step)




    contains

        ! allocate memory of temporaries
        subroutine classic_rk4_alloc(num)
            integer, intent(in) :: num

            allocate(inipos(num, 2))
            allocate(inivor(num, 1))

            allocate(k1o(num, 2))
            allocate(k2o(num, 2))
            allocate(k3o(num, 2))
            allocate(k4o(num, 2))

            allocate(w1o(num, 1))
            allocate(w2o(num, 1))
            allocate(w3o(num, 1))
            allocate(w4o(num, 1))

            allocate(strain(num, 4))

            if (parcel_info%is_elliptic) then
                allocate(iniB(num, 2))

                allocate(b1o(num, 2))
                allocate(b2o(num, 2))
                allocate(b3o(num, 2))
                allocate(b4o(num, 2))
            endif

        end subroutine classic_rk4_alloc

        ! deallocate memory of temporaries
        subroutine classic_rk4_dealloc

            deallocate(inipos)
            deallocate(inivor)

            deallocate(k1o)
            deallocate(k2o)
            deallocate(k3o)
            deallocate(k4o)

            deallocate(w1o)
            deallocate(w2o)
            deallocate(w3o)
            deallocate(w4o)

            deallocate(strain)

            if (parcel_info%is_elliptic) then
                deallocate(iniB)

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

            ! copy input position and B matrix
            inipos(1:n_parcels,:) = parcels%position(1:n_parcels,:)
            inivor(1:n_parcels,:) = parcels%vorticity(1:n_parcels,:)
            iniB(1:n_parcels,:) = parcels%B(1:n_parcels,:)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k1o, w1o, strain)
            b1o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k1o)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) + 0.5_dp * dt * k1o(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + 0.5_dp * dt * w1o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + 0.5_dp * dt * b1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k1o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k2o, w2o, strain)
            b2o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k2o)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) + 0.5_dp * dt * k2o(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + 0.5_dp * dt * w2o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + 0.5_dp * dt * b2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k2o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k3o, w3o, strain)
            b3o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))


            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k3o)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) + dt * k3o(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + dt * w3o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + dt * b3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k3o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k4o, w4o, strain)
            b4o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k4o)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                             + dt / 6.0_dp * (k1o(1:n_parcels,:) + 2.0_dp &
                             * k2o(1:n_parcels,:) + 2.0_dp * k3o(1:n_parcels,:) + k4o(1:n_parcels,:))

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                             + dt / 6.0_dp * (w1o(1:n_parcels,:) + 2.0_dp &
                             * w2o(1:n_parcels,:) + 2.0_dp * w3o(1:n_parcels,:) + w4o(1:n_parcels,:))

            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) &
                      + dt / 6.0_dp * (b1o(1:n_parcels,:) + 2.0_dp * b2o(1:n_parcels,:) &
                      + 2.0_dp * b3o(1:n_parcels,:) + b4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

        end subroutine classic_rk4_elliptic


        subroutine classic_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            ! copy input position and B matrix
            inipos(1:n_parcels,:) = parcels%position(1:n_parcels,:)
            inivor(1:n_parcels,:) = parcels%vorticity(1:n_parcels,:)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k1o, w1o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k1o)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                          + 0.5_dp * dt * k1o(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                          + 0.5_dp * dt * w1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k1o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k2o, w2o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k2o)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                          + 0.5_dp * dt * k2o(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                          + 0.5_dp * dt * w2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k2o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k3o, w3o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k3o)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + dt * k3o(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                            + dt * w3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, k3o)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call grid2par(k4o, w4o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, k4o)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                             + dt / 6.0_dp * (k1o(1:n_parcels,:) + 2.0_dp &
                             * k2o(1:n_parcels,:) + 2.0_dp * k3o(1:n_parcels,:) + k4o(1:n_parcels,:))

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                             + dt / 6.0_dp * (w1o(1:n_parcels,:) + 2.0_dp &
                             * w2o(1:n_parcels,:) + 2.0_dp * w3o(1:n_parcels,:) + w4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, k4o)

        end subroutine classic_rk4_non_elliptic

end module classic_rk4
