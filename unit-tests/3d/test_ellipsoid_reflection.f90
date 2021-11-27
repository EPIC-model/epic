! =============================================================================
!                       Test ellipsoid reflection
!
!         This unit test checks the reflection of parcels which is
!         applied after splitting to ensure the centres of the child
!         parcels are within the domain.
! =============================================================================
program test_ellipsoid_reflection
    use unit_test
    use options, only : verbose
    use constants, only : pi, zero, one, two, three, six, f12, f13, f14, f18, f23, twopi
    use parcel_container, only : parcel_alloc   &
                               , parcel_dealloc &
                               , n_parcels      &
                               , parcels
    use parcel_bc, only : apply_reflective_bc
    use parcel_ellipsoid, only : get_abc        &
                               , get_angles     &
                               , get_B33
    use parameters, only : update_parameters, lower, dx, vcell, extent, nx, ny, nz
    implicit none

    integer :: iter_t, iter_p
    double precision :: abc, abc23, B11, B12, B13, B22, B23, error, a2, b2, c2
    double precision :: st, sp, ct, cp, theta, phi, angle

    call parse_command_line

    nx = 2
    ny = 2
    nz = 2
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

    if (verbose) then
        open(80, file='original_position_and_Bmatrix.asc')
        open(81, file='reflected_position_and_Bmatrix.asc')
    endif

    do iter_t = 0, 359
        do iter_p = 0, 179
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

            if (verbose) then
                write(80, *) parcels%position(1, :), parcels%B(1, :), &
                             get_B33(parcels%B(1, :), parcels%volume(1))
            endif

            call apply_reflective_bc(parcels%position(1, :), parcels%B(1, :))

            if (verbose) then
                write(81, *) parcels%position(1, :), parcels%B(1, :), &
                             get_B33(parcels%B(1, :), parcels%volume(1))
            endif

            call check_result('lower')

            !
            ! upper boundary check
            !
            parcels%position(1, :) = (/f12, f12, one + f12 * dx(3)/)
            parcels%B(1, 1) = B11
            parcels%B(1, 2) = B12
            parcels%B(1, 3) = B13
            parcels%B(1, 4) = B22
            parcels%B(1, 5) = B23

            call apply_reflective_bc(parcels%position(1, :), parcels%B(1, :))

            call check_result('upper')
        enddo
    enddo

    if (verbose) then
        close(80)
        close(81)
    endif

    call print_result_dp('Test ellipsoid reflection', error, atol=3.0e-13)

    call parcel_dealloc


    contains

        subroutine check_result(bc)
            character(*), intent(in) :: bc
            double precision         :: angles(2)

            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(1, 2) - B12))
            error = max(error, abs(parcels%B(1, 3) + B13))
            error = max(error, abs(parcels%B(1, 4) - B22))
            error = max(error, abs(parcels%B(1, 5) + B23))

            ! (/azimuth, polar/)
            angles = get_angles(parcels%B(1, :), parcels%volume(1))

            ! -pi <= theta <= pi
            if (angles(1) < 0) then
                angles(1) = twopi + angles(1)
            endif

            ! rotated by pi/2
            if (dabs(dabs(theta - angles(1)) - pi) < 10000.0d0 * epsilon(pi)) then
                if (theta > angles(1)) then
                    angles(1) = pi + angles(1)
                else
                    angles(1) = angles(1) - pi
                endif
            endif


            ! 0 <= angles(2) <= pi/2
            if (phi > pi/2) then
                angles(2) = pi - dabs(angles(2))
            endif

            ! Neglect if phi = 0, pi/2 and pi, since then theta does not matter
            if ((iter_p .ne. 0) .and. (iter_p .ne. 90) .and. (iter_p .ne. 180)) then
                error = max(error, abs(theta - angles(1)))
            endif

            error = max(error, dabs(dabs(phi) - dabs(angles(2))))

            error = max(error, abs(f12 - parcels%position(1, 1)))
            error = max(error, abs(f12 - parcels%position(1, 2)))

            if (bc == 'lower') then
                error = max(error, abs(f12 * dx(3) - parcels%position(1, 3)))
            else
                error = max(error, abs(one - f12 * dx(3) - parcels%position(1, 3)))
            endif
        end subroutine check_result


end program test_ellipsoid_reflection
