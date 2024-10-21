! =============================================================================
!                       Test ellipse reflection
!
!         This unit test checks the reflection of parcels which is
!         applied after splitting to ensure the centres of the child
!         parcels are within the domain.
! =============================================================================
program test_ellipse_reflection
    use unit_test
    use constants, only : pi, zero, one, two, three, four, f12, f14
    use parcel_container
    use parcel_bc, only : apply_reflective_bc
    use parcel_ellipse, only : get_ab, get_angle, get_B22
    use parameters, only : update_parameters, lower, dx, vcell, extent, nx, nz
    implicit none

    integer :: iter
    double precision :: angle, ab, B11, B12, B22, error, a2, b2

    nx = 10
    nz = 10
    extent = (/one, one/)
    lower  = (/zero, zero/)

    call update_parameters

    n_parcels = 1
    call parcel_alloc(1)

    error = zero

    parcels%volume(1) = f14 * vcell

    ab = get_ab(parcels%volume(1))
    a2 = four * ab
    b2 = f14 * ab

    do iter = 0, 360
        angle = dble(iter) * pi / 180.0d0

        !
        ! lower boundary check
        !

        parcels%position(:, 1) = (/f12, - f12 * dx(2)/)
        B11 = a2 * dcos(angle) ** 2 + b2 * dsin(angle) ** 2
        B12 = f12 * (a2 - b2) * dsin(two * angle)
        parcels%B(1, 1) = B11
        parcels%B(2, 1) = B12

        B22 = get_B22(B11, B12, V)
        parcels%B(3, 1) = B22

        call apply_reflective_bc(parcels%position(:, 1), parcels%B(:, 1))

        call check_result('lower')

        !
        ! upper boundary check
        !

        angle = dble(iter) * pi / 180.0d0

        parcels%position(:, 1) = (/f12, one + f12 * dx(2)/)
        B11 = a2 * dcos(angle) ** 2 + b2 * dsin(angle) ** 2
        B12 = f12 * (a2 - b2) * dsin(two * angle)
        parcels%B(1, 1) = B11
        parcels%B(2, 1) = B12

        call apply_reflective_bc(parcels%position(:, 1), parcels%B(:, 1))

        call check_result('upper')
    enddo

    call print_result_dp('Test ellipse reflection', error, atol=1.0e-14)

    call parcel_dealloc


    contains

        subroutine check_result(bc)
            character(*), intent(in) :: bc

            ! get_angle computes the angle in the first and fourth quadrant, i.e.,
            ! -pi/2 <= get_angle <= pi/2
            if (angle > pi / two .and. angle <= three * pi / two) then
                angle = angle - pi
            else if (angle > three * pi / two) then
                angle = angle - two * pi
            endif

            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(2, 1) + B12))
            error = max(error, abs(angle + get_angle(parcels%B(1, 1), parcels%B(2, 1), &
                                                     parcels%volume(1))))
            error = max(error, abs(f12 - parcels%position(1, 1)))

            if (bc == 'lower') then
                error = max(error, abs(f12 * dx(2) - parcels%position(2, 1)))
            else
                error = max(error, abs(one - f12 * dx(2) - parcels%position(2, 1)))
            endif
        end subroutine check_result


end program test_ellipse_reflection
