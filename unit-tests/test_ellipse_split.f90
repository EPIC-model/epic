! =============================================================================
!                       Test ellipse split
!
!         This unit test checks the splitting of an ellipse. The ellipse
!         is centred at the origin with an orientation of 45 degrees.
! =============================================================================
program test_ellipse_split
    use unit_test
    use constants, only : pi, zero, one, three, four, five, ten
    use parcel_container
    use parcel_split, only : split_ellipses
    use parameters, only : update_parameters, nx, nz, extent, lower
    implicit none

    double precision, parameter :: lam = five
    double precision, parameter :: angle = 0.25d0 * pi
    double precision, parameter :: evec(2) = (/dcos(angle), dsin(angle)/)
    double precision :: h, ab, B11, B12, B22, pos(2, 2), error, a2, b2

    nx = 10
    nz = 10
    extent = (/ten, ten/)
    lower = (/-five, -five/)
    call update_parameters

    n_parcels = 1
    call parcel_alloc(2)

    ab = one

    a2 = ab * lam
    b2 = ab / lam

    parcels%position(1, :) = zero
    parcels%volume(1) = ab * pi
    parcels%buoyancy(1) = one
    parcels%humidity(1) = one

    B11 = a2 * dcos(angle) ** 2 + b2 * dsin(angle) ** 2
    B12 = 0.5d0 * (a2 - b2) * dsin(2.0 * angle)
    B22 = a2 * dsin(angle) ** 2 + b2 * dcos(angle) ** 2

    parcels%B(1, 1) = B11
    parcels%B(1, 2) = B12

    ! analytic split
    h = 0.25d0 * dsqrt(three * a2)
    B11 = B11 - 0.75d0 * a2 * evec(1) ** 2
    B12 = B12 - 0.75d0 * a2 * evec(1) * evec(2)
    pos(1, :) = parcels%position(1, :) + h * evec
    pos(2, :) = parcels%position(1, :) - h * evec

    ! numerical split
    call split_ellipses(parcels, threshold=four, vthreshold=one)

    !
    ! check result
    !

    error = zero

    ! first parcel
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, sum(abs(pos(1, :) - parcels%position(1, :))))
    error = max(error, abs(0.5d0 * ab * pi - parcels%volume(1)))
    error = max(error, abs(parcels%buoyancy(1) - one))
    error = max(error, abs(parcels%humidity(1) - one))


    ! second parcel
    error = max(error, abs(parcels%B(2, 1) - B11))
    error = max(error, abs(parcels%B(2, 2) - B12))
    error = max(error, sum(abs(pos(2, :) - parcels%position(2, :))))
    error = max(error, abs(0.5d0 * ab * pi - parcels%volume(2)))
    error = max(error, dble(abs(n_parcels - 2)))
    error = max(error, abs(parcels%buoyancy(2) - one))
    error = max(error, abs(parcels%humidity(2) - one))

    call print_result_dp('Test ellipse split', error)

    call parcel_dealloc

end program test_ellipse_split
