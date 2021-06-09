! =============================================================================
!               Classical 4th order Runge-Kutta method
!     (see https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
! =============================================================================
module classic_rk4
    use constants, only : max_num_parcels, f16
    use parameters, only : nx, nz
    use options, only : parcel
    use parcel_container
    use parcel_bc
    use rk4_utils, only: get_B
    use tri_inversion, only : vor2vel, vorticity_tendency
    use parcel_interpl, only : par2grid, grid2par
    use fields, only : velgradg, velog, vortg, vtend, tbuoyg
    implicit none
    integer, parameter :: dp=kind(0.d0)           ! double precision

    ! classic_rk4 temporaries
    double precision, allocatable, dimension(:, :) :: &
        veo,                    &   ! sum of velocities
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

            allocate(veo(num, 2))

            allocate(w1o(num, 1))
            allocate(w2o(num, 1))
            allocate(w3o(num, 1))
            allocate(w4o(num, 1))

            allocate(strain(num, 4))

            if (parcel%is_elliptic) then
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

            deallocate(veo)

            deallocate(w1o)
            deallocate(w2o)
            deallocate(w3o)
            deallocate(w4o)

            deallocate(strain)

            if (parcel%is_elliptic) then
                deallocate(iniB)

                deallocate(b1o)
                deallocate(b2o)
                deallocate(b3o)
                deallocate(b4o)
            endif

        end subroutine classic_rk4_dealloc


        subroutine classic_rk4_step(dt)
            double precision, intent(in) :: dt

            if (parcel%is_elliptic) then
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
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w1o, strain)
            b1o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + 0.5_dp * dt * parcels%velocity(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + 0.5_dp * dt * w1o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + 0.5_dp * dt * b1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w2o, strain)
            b2o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + 0.5_dp * dt * parcels%velocity(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + 0.5_dp * dt * w2o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + 0.5_dp * dt * b2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = veo(1:n_parcels, :) + 2.0_dp * parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w3o, strain)
            b3o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))


            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + dt * parcels%velocity(1:n_parcels,:)
            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) + dt * w3o(1:n_parcels,:)
            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) + dt * b3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = veo(1:n_parcels, :) + 2.0_dp * parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w4o, strain)
            b4o(1:n_parcels,:) = get_B(parcels%B(1:n_parcels,:), strain(1:n_parcels,:), &
                                       parcels%volume(1:n_parcels))

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)


            parcels%velocity(1:n_parcels, :) = (veo(1:n_parcels, :) + parcels%velocity(1:n_parcels, :)) * f16

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + dt * parcels%velocity(1:n_parcels, :)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                             + dt * f16 * (w1o(1:n_parcels,:) + 2.0_dp &
                             * w2o(1:n_parcels,:) + 2.0_dp * w3o(1:n_parcels,:) + w4o(1:n_parcels,:))

            parcels%B(1:n_parcels,:) = iniB(1:n_parcels, :) &
                      + dt * f16 * (b1o(1:n_parcels,:) + 2.0_dp * b2o(1:n_parcels,:) &
                      + 2.0_dp * b3o(1:n_parcels,:) + b4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_elliptic


        subroutine classic_rk4_non_elliptic(dt)
            double precision, intent(in) :: dt

            ! copy input position and B matrix
            inipos(1:n_parcels,:) = parcels%position(1:n_parcels,:)
            inivor(1:n_parcels,:) = parcels%vorticity(1:n_parcels,:)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w1o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)


            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                          + 0.5_dp * dt * parcels%velocity(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                          + 0.5_dp * dt * w1o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w2o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                          + 0.5_dp * dt * parcels%velocity(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                          + 0.5_dp * dt * w2o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = veo(1:n_parcels, :) &
                                + 2.0_dp * parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w3o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                                            + dt * parcels%velocity(1:n_parcels,:)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                                            + dt * w3o(1:n_parcels,:)

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)
            veo(1:n_parcels, :) = veo(1:n_parcels, :) &
                                + 2.0_dp * parcels%velocity(1:n_parcels, :)

            call par2grid
            call vor2vel(vortg, velog, velgradg)
            call vorticity_tendency(tbuoyg, vtend)
            call grid2par(parcels%velocity, w4o, strain)

            ! apply velocity BC --> only important for free slip
            call apply_parcel_bc(parcels%position, parcels%velocity)

            parcels%velocity(1:n_parcels, :) = (veo(1:n_parcels, :) + parcels%velocity(1:n_parcels, :)) * f16

            parcels%position(1:n_parcels,:) = inipos(1:n_parcels, :) &
                             + dt * parcels%velocity(1:n_parcels, :)

            parcels%vorticity(1:n_parcels,:) = inivor(1:n_parcels, :) &
                             + dt * f16 * (w1o(1:n_parcels,:) + 2.0_dp &
                             * w2o(1:n_parcels,:) + 2.0_dp * w3o(1:n_parcels,:) + w4o(1:n_parcels,:))

            ! apply position BC
            call apply_parcel_bc(parcels%position, parcels%velocity)

        end subroutine classic_rk4_non_elliptic

end module classic_rk4
