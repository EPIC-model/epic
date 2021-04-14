! =============================================================================
!                       Test ellipse orientation
!
!       This unit test checks the ellipse orientation computed using the
!       B matrix values, i.e., B11 and B12. The angle computed via the
!       B matrix is within [-pi/2, pi/2].
! =============================================================================
program ellipse_orientation
    use constants, only : pi
    use parcel_container, only : parcels, n_parcels, parcel_alloc, parcel_dealloc
    use ellipse
    implicit none

    double precision :: origin(2) = (/0.0, 0.0/)
    double precision :: extent(2) =  (/0.2, 0.2/)
    integer :: iter
    integer :: grid(2) = (/2, 2/)
    double precision :: angle, B11, B12, V
    double precision, parameter :: lam = 3.0
    logical :: failed = .false.

    n_parcels = 1
    call parcel_alloc(1)

    parcels%position = 0.0
    parcels%volume = 0.25 * product(extent / (grid - 1))

    V = parcels%volume(1, 1)

    do iter = 0, 360

        angle = dble(iter) * pi / 180.0d0

        B11 = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2

        B12 = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)

        parcels%B(1, 1) = B11

        parcels%B(1, 2) = B12

        ! get_angle computes the angle in the first and fourth quadrant, i.e.,
        ! -pi/2 <= get_angle <= pi/2
        if (angle > pi / 2 .and. angle <= 3.0 * pi / 2) then
            angle = angle - pi
        else if (angle > 3.0 * pi / 2) then
            angle = angle - 2 * pi
        endif

        failed = (failed .or. abs(angle - get_angle(B11, B12, V)) > 1.0e-13)

    enddo

    call parcel_dealloc()

    if (failed) then
        print '(a26, a10)', 'Test ellipse orientation:', 'FAILED'
    else
        print '(a26, a10)', 'Test ellipse orientation:', 'PASSED'
    endif

end program ellipse_orientation
