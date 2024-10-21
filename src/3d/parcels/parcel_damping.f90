! =============================================================================
! This module contains the subroutines to do damping of parcel properties to
! gridded fields in a conservative manner. It does this by nudging all the
! parcels associated with a grid point to the grid point value, with the strength
! of damping proportional to the grid point contribution to the gridded value and
! the strain rate at the grid point.
! =============================================================================
module parcel_damping
    use constants, only :  f14, zero, one
    use mpi_timer, only : start_timer, stop_timer
    use parameters, only : nx, nz, vmin
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid
    use parcel_interpl
    use fields
    use parameters, only : lower, upper
    use omp_lib
    use mpi_layout, only : box
    use mpi_utils, only : mpi_exit_on_error
    use options, only : damping
    use inversion_mod, only : vor2vel
    use rk_utils, only : get_strain_magnitude_field
    implicit none


    private

    ! interpolation indices
    ! (first dimension x, y, z; second dimension l-th index)
    integer :: is, js, ks
    integer :: damping_timer

    ! interpolation weights
    double precision :: weights(0:1,0:1,0:1)
    double precision :: weight(0:1,0:1,0:1)
    double precision :: time_fact(0:1,0:1,0:1)

    public :: parcel_damp, damping_timer

    contains

        subroutine parcel_damp(dt)
            double precision, intent(in)  :: dt

            if (damping%l_vorticity .and. damping%l_surface_vorticity) then
                call mpi_exit_on_error("damping%l_vorticity and  damping%l_surface_vorticity both activated, only one allowed")
            elseif (damping%l_scalars .and. damping%l_surface_scalars) then
                call mpi_exit_on_error("damping%l_scalars and  damping%l_surface_scalars both activated, only one allowed")
            endif
            if (damping%l_vorticity .or. damping%l_scalars .or. &
                damping%l_surface_vorticity .or. damping%l_surface_scalars) then
                ! ensure gridded fields are up to date
                call par2grid(.false.)
                call vor2vel

                ! Reflect beyond boundaries to ensure damping is conservative
                ! This is because the points below the surface contribute to the level above
                !$omp parallel workshare
                vortg(-1,   :, :, :) = vortg(1, :, :, :)
                vortg(nz+1, :, :, :) = vortg(nz-1, :, :, :)

#ifndef ENABLE_DRY_MODE
                humg(-1,   :, :) = humg(1, :, :)
                humg(nz+1, :, :) = humg(nz-1, :, :)
                dbuoyg(-1,   :, :) = dbuoyg(1, :, :)
                dbuoyg(nz+1, :, :) = dbuoyg(nz-1, :, :)
#else
                tbuoyg(-1,   :, :) = tbuoyg(1, :, :)
                tbuoyg(nz+1, :, :) = tbuoyg(nz-1, :, :)
#endif
                !$omp end parallel workshare

                call get_strain_magnitude_field
                call perturbation_damping(dt, .true.)
            end if

        end subroutine parcel_damp

        !
        ! @pre: the strain must be calculated and the gridded fields updated
        subroutine perturbation_damping(dt, l_reuse)
            double precision, intent(in)  :: dt
            logical, intent(in)           :: l_reuse
#if defined (ENABLE_P2G_1POINT) && !defined (NDEBUG)
            logical                       :: l_reuse_dummy
#endif
            integer                       :: n, p, l, surface_index
            double precision              :: points(3, n_points_p2g)
            double precision              :: pvol
            ! tendencies need to be summed up between associated 4 points
            ! before modifying the parcel attribute
            double precision              :: vortend(3)
#ifndef ENABLE_DRY_MODE
            double precision              :: dbuoytend
            double precision              :: humtend
#else
            double precision              :: tbuoytend
#endif
            call start_timer(damping_timer)

            ! This is only here to allow debug compilation
            ! with a warning for unused variables
#if defined (ENABLE_P2G_1POINT) && !defined (NDEBUG)
            l_reuse_dummy=l_reuse
#endif

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, weight, surface_index) &
#ifndef ENABLE_DRY_MODE
            !$omp& private(is, js, ks, weights, vortend, humtend, dbuoytend, time_fact)
#else
            !$omp& private(is, js, ks, weights, vortend, tbuoytend, time_fact)
#endif
            do n = 1, n_parcels
                ! check if only surface damping applies and we are far from surfaces
                ! put in a buffer here as parcels can get stretched in integration
                if(.not.(damping%l_vorticity .or. damping%l_scalars)) then
                    if(parcels%position(3, n) > lower(3) + 2 * dx(3)) then
                    if(parcels%position(3, n) < upper(3) - 2 * dx(3)) then
                        cycle 
                    end if
                    end if
                endif

                pvol = parcels%volume(n)
#ifndef ENABLE_P2G_1POINT
                points = get_ellipsoid_points(parcels%position(:, n), &
                                              parcels%B(:, n), n, l_reuse)
#else
                points(:, 1) = parcels%position(:, n)
#endif
                vortend = zero
#ifndef ENABLE_DRY_MODE
                humtend = zero
                dbuoytend = zero
#else
                tbuoytend = zero
#endif

                ! we have 4 points per ellipsoid
                do p = 1, n_points_p2g
                    call trilinear(points(:, p), is, js, ks, weights)
                    weight = point_weight_p2g * weights

                    if (damping%l_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
                        do l = 1, 3
                            vortend(l) = vortend(l)+sum(weight * time_fact * (vortg(ks:ks+1, js:js+1, is:is+1, l) &
                                       - parcels%vorticity(l,n)))
                        enddo
                    endif

                    if (damping%l_scalars) then
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
#ifndef ENABLE_DRY_MODE
                        humtend = humtend + sum(weight * time_fact * (humg(ks:ks+1, js:js+1, is:is+1) - parcels%humidity(n)))
                        dbuoytend = dbuoytend + sum(weight * time_fact * (dbuoyg(ks:ks+1, js:js+1, is:is+1) - parcels%buoyancy(n)))
#else
                        tbuoytend = tbuoytend + sum(weight * time_fact * (tbuoyg(ks:ks+1, js:js+1, is:is+1) - parcels%buoyancy(n)))
#endif
                    endif
                    
                    if (damping%l_surface_vorticity .or. damping%l_surface_scalars) then
                        ! Index to keep track of grid cells right above/below boundary
                        ! This is because the damping only happens at the boundary level
                        ! Consistent with reflection used in parcel_damp
                        if ((ks == box%lo(3)-1) .or. (ks == box%hi(3)-1)) then
                            surface_index = 1 ! below lower or below upper boundary
                        elseif ((ks == box%lo(3)) .or. (ks == box%hi(3))) then
                            surface_index = 0 ! above lower or above upper boundary
                        else
                            cycle ! continue loop if not near a surface
                        endif
                    else
                        cycle ! continue loop if no surface damping
                    endif    

                    if (damping%l_surface_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
                        do l = 1, 3
                            vortend(l) = vortend(l)+sum(weight(surface_index, :, :) * time_fact(surface_index, :, :) * &
                                         (vortg(ks+surface_index, js:js+1, is:is+1, l)  - parcels%vorticity(l,n)))
                        enddo
                    endif

                    if (damping%l_surface_scalars) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
#ifndef ENABLE_DRY_MODE
                        humtend = humtend + sum(weight(surface_index, :, :) * time_fact(surface_index, :, :) * &
                                                  (humg(ks+surface_index, js:js+1, is:is+1) - parcels%humidity(n)))
                        dbuoytend = dbuoytend + sum(weight(surface_index, :, :) * time_fact(surface_index, :, :) * &
                                                  (dbuoyg(ks+surface_index, js:js+1, is:is+1) - parcels%buoyancy(n)))
#else
                        tbuoytend = tbuoytend + sum(weight(surface_index, :, :) * time_fact(surface_index, :, :) * &
                                                        (tbuoyg(ks+surface_index, js:js+1, is:is+1) - parcels%buoyancy(n)))
#endif
                    endif
                enddo
                ! Add all the tendencies only at the end
                if (damping%l_vorticity .or. damping%l_surface_vorticity) then
                    do l=1, 3
                        parcels%vorticity(l,n) = parcels%vorticity(l,n) + vortend(l)
                    enddo
                endif
                if (damping%l_scalars .or. damping%l_surface_scalars) then
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n) = parcels%humidity(n) + humtend
                    parcels%buoyancy(n) = parcels%buoyancy(n) + dbuoytend
#else
                    parcels%buoyancy(n) = parcels%buoyancy(n) + tbuoytend
#endif
                endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(damping_timer)

        end subroutine perturbation_damping


end module parcel_damping
