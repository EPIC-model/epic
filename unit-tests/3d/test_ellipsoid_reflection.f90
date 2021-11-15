! =============================================================================
!                       Test ellipsoid reflection
!
!         This unit test checks the reflection of parcels which is
!         applied after splitting to ensure the centres of the child
!         parcels are within the domain.
! =============================================================================
program test_ellipsoid_reflection
    use unit_test
    use constants, only : pi, zero, one, two, three, six, f12, f13, f18, f23
    use parcel_container, only : parcel_alloc   &
                               , parcel_dealloc &
                               , n_parcels      &
                               , parcels
    use parcel_bc, only : apply_reflective_bc
    use parcel_ellipsoid, only : get_abc             &
                               , get_azimuthal_angle &
                               , get_polar_angle
    use parameters, only : update_parameters, lower, dx, vcell, extent, nx, ny, nz
    implicit none

    integer :: iter_t, iter_p
    double precision :: abc, abc23, B11, B12, B13, B22, B23, error, a2, b2, c2
    double precision :: st, sp, ct, cp, theta, phi, angle

    nx = 10
    ny = 10
    nz = 10
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)

    call update_parameters

    n_parcels = 1
    call parcel_alloc(1)

    error = zero

    parcels%volume(1) = f18 * vcell

    abc = get_abc(parcels%volume(1))
    abc23 = abc ** f23
    a2 = six * abc23
    b2 = f12 * abc23
    c2 = f13 * abc23

    do iter_t = 0, 360
        do iter_p = 0, 180
            theta = dble(iter_t) * pi / 180.0d0     ! azimuthal angle, [0, 2pi[
            phi = dble(iter_p) * pi / 180.0d0       ! polar angle, [0, pi]

            st = dsin(theta)
            ct = dcos(theta)

            sp = dsin(phi)
            cp = dcos(phi)

            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp
            ! B33 = a2 * cp ** 2 + c2 * sp ** 2

            !
            ! lower boundary check
            !
            parcels%position(1, :) = (/f12, f12, - f12 * dx(3)/)
            parcels%B(1, 1) = B11
            parcels%B(1, 2) = B12
            parcels%B(1, 3) = B13
            parcels%B(1, 4) = B22
            parcels%B(1, 5) = B23

            angle = get_azimuthal_angle(parcels%B(1, :), parcels%volume(1))
            print *, "before", angle
!             stop

            call apply_reflective_bc(parcels%position(1, :), parcels%B(1, :))

            call check_result('lower')

!             !
!             ! upper boundary check
!             !
!             parcels%position(1, :) = (/f12, f12, one + f12 * dx(3)/)
!             parcels%B(1, 1) = B11
!             parcels%B(1, 2) = B12
!             parcels%B(1, 3) = B13
!             parcels%B(1, 4) = B22
!             parcels%B(1, 5) = B23
!
!             call apply_reflective_bc(parcels%position(1, :), parcels%B(1, :))
!
!             call check_result('upper')
        enddo
    enddo

    call print_result_dp('Test ellipsoid reflection', error, atol=1.0e-14)

    call parcel_dealloc


    contains

        subroutine check_result(bc)
            character(*), intent(in) :: bc
            double precision         :: azimuth, polar

            ! get_angle computes the angle in the first and fourth quadrant, i.e.,
            ! -pi/2 <= get_angle <= pi/2
!             if (theta > pi / two .and. theta <= three * pi / two) then
!                 theta = theta - pi
!             else if (theta > three * pi / two) then
!                 theta = pi - theta
!             endif
            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(1, 2) - B12))
            error = max(error, abs(parcels%B(1, 3) + B13))
            error = max(error, abs(parcels%B(1, 4) - B22))
            error = max(error, abs(parcels%B(1, 5) + B23))

            azimuth = get_azimuthal_angle(parcels%B(1, :), parcels%volume(1))
            polar = get_polar_angle(parcels%B(1, :), parcels%volume(1))

            if (azimuth < zero) then
                azimuth = pi + azimuth
            endif

            if (dabs(dabs(theta - azimuth) - pi) <= 10.0d0 * epsilon(pi)) then
                azimuth = pi - azimuth
                print *, "HI"
            endif

            error = max(error, abs(theta - azimuth))

            print *, theta, azimuth, phi
            if (error > 1.0e-14) then
                print *, "azimuth", azimuth
                print *, "polar", polar
                print *, "theta", theta
                print *, "phi", phi
                print *, "B11", B11
                print *, "B12", B12
                print *, "B13", B13
                print *, "B22", B22
                print *, "B23", B23
                print *, "B33", a2 * cp ** 2 + c2 * sp ** 2
                print *, "a2", a2
                print *, "b2", b2
                print *, "c2", c2
                print *, "st", st
                print *, "ct", ct
                print *, "sp", sp
                print *, "cp", cp
                stop
            endif
!             error = max(error, abs(phi + get_polar_angle(parcels%B(1, :), &
!                                                          parcels%volume(1))))
            error = max(error, abs(f12 - parcels%position(1, 1)))
            error = max(error, abs(f12 - parcels%position(1, 2)))

            if (bc == 'lower') then
                error = max(error, abs(f12 * dx(3) - parcels%position(1, 3)))
            else
                error = max(error, abs(one - f12 * dx(3) - parcels%position(1, 3)))
            endif
        end subroutine check_result


end program test_ellipsoid_reflection
