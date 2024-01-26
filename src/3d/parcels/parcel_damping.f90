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
    use parcels_mod, only : parcels
    use parcel_interpl
    use fields
    use omp_lib
    use options, only : damping
    use interpl, only : trilinear
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

            if (damping%l_vorticity .or. damping%l_scalars) then
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
#endif
                tbuoyg(-1,   :, :) = tbuoyg(1, :, :)
                tbuoyg(nz+1, :, :) = tbuoyg(nz-1, :, :)
                !$omp end parallel workshare

                call get_strain_magnitude_field
                call perturbation_damping(dt, .true.)
            end if

        end subroutine parcel_damp

        !
        ! @pre
        subroutine perturbation_damping(dt, l_reuse)
            double precision, intent(in)  :: dt
            logical, intent(in)           :: l_reuse
            integer                       :: n, p, l
            double precision              :: points(3, n_points_p2g)
            double precision              :: pvol
            ! tendencies need to be summed up between associated 4 points
            ! before modifying the parcel attribute
            double precision              :: vortend(3)
            double precision              :: buoytend
#ifndef ENABLE_DRY_MODE
            double precision              :: humtend
#endif

            call start_timer(damping_timer)

            !$omp parallel default(shared)
            !$omp do private(n, p, l, points, pvol, weight) &
#ifndef ENABLE_DRY_MODE
            !$omp& private(is, js, ks, weights, vortend, buoytend, humtend, time_fact)
#else
            !$omp& private(is, js, ks, weights, vortend, buoytend, time_fact)
#endif
            do n = 1, parcels%local_num
                pvol = parcels%volume(n)
#ifndef ENABLE_P2G_1POINT
                points = parcels%get_points(n, l_reuse)
#else
                points(:, 1) = parcels%position(:, n)
#endif
                vortend = zero
                buoytend = zero
#ifndef ENABLE_DRY_MODE
                humtend = zero
#endif

                ! we have 4 points per ellipsoid
                do p = 1, n_points_p2g
                    call trilinear(points(:, p), is, js, ks, weights)
                    weight = point_weight_p2g * weights
                    if (damping%l_vorticity) then
                        ! Note this exponential factor can be different for vorticity/scalars
                        time_fact = one - exp(-damping%vorticity_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
                        do l = 1,3
                            vortend(l) = vortend(l)+sum(weight * time_fact * (vortg(ks:ks+1, js:js+1, is:is+1, l) &
                                       - parcels%vorticity(l,n)))
                        enddo
                    endif
                    if (damping%l_scalars) then
                        time_fact = one - exp(-damping%scalars_prefactor * strain_mag(ks:ks+1, js:js+1, is:is+1) * dt)
#ifndef ENABLE_DRY_MODE
                        humtend = humtend + sum(weight * time_fact * (humg(ks:ks+1, js:js+1, is:is+1) - parcels%humidity(n)))
                        buoytend = buoytend + sum(weight * time_fact * (dbuoyg(ks:ks+1, js:js+1, is:is+1) - parcels%buoyancy(n)))
#else
                        buoytend = buoytend + sum(weight * time_fact * (tbuoyg(ks:ks+1, js:js+1, is:is+1) - parcels%buoyancy(n)))
#endif
                    endif
                enddo
                ! Add all the tendencies only at the end
                if (damping%l_vorticity) then
                    do l=1,3
                        parcels%vorticity(l,n) = parcels%vorticity(l,n) + vortend(l)
                    enddo
                endif
                if (damping%l_scalars) then
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n) = parcels%humidity(n) + humtend
#endif
                    parcels%buoyancy(n) = parcels%buoyancy(n) + buoytend
                endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(damping_timer)

        end subroutine perturbation_damping


end module parcel_damping
