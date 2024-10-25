! =============================================================================
!                       Test ellipse orientation
!
!       This unit test checks the ellipse orientation computed using the
!       B matrix values, i.e., B11 and B12. The angle computed via the
!       B matrix is within [-pi/2, pi/2].
! =============================================================================
program test_ellipse_orientation
    use unit_test
    use constants, only : pi, zero, two, three, f14
    use parcel_container, only : parcels, n_parcels, parcel_alloc, parcel_dealloc
    use parcel_ellipse
    implicit none

    double precision :: extent(2) =  (/0.2d0, 0.2d0/)
    integer :: iter
    integer :: grid(2) = (/2, 2/)
    double precision :: angle, B11, B12, V, a2, b2
    double precision, parameter :: lam = three
    logical :: passed = .true.

    n_parcels = 1
    call parcel_alloc(1)

    parcels%position = zero
    parcels%volume = f14 * product(extent / (grid - 1))

    V = parcels%volume(1)

    a2 = (lam * V) / pi
    b2 = V / (pi * lam)

    do iter = 0, 360

        angle = dble(iter) * pi / 180.0d0

        B11 = a2 * cos(angle) ** 2 + b2 * sin(angle) ** 2

        B12 = (a2 - b2) * sin(angle) * cos(angle)

        parcels%B(1, 1) = B11

        parcels%B(2, 1) = B12

        ! get_angle computes the angle in the first and fourth quadrant, i.e.,
        ! -pi/2 <= get_angle <= pi/2
        if (angle > pi / two .and. angle <= three * pi / two) then
            angle = angle - pi
        else if (angle > three * pi / two) then
            angle = angle - two * pi
        endif

        passed = (passed .and. abs(angle - get_angle(B11, B12, V)) < dble(1.0e-13))

    enddo

    call parcel_dealloc()

    call print_result_logical('Test ellipse orientation', passed)

end program test_ellipse_orientation
