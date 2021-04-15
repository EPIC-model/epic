! =============================================================================
!                       Test ellipse split
!
!         This unit test checks the splitting of an ellipse. The ellipse
!         is centred at the origin with an orientation of 45 degrees.
! =============================================================================
program test_ellipse_split
    use constants, only : pi
    use parcel_container
    use parcel_split, only : split_ellipses
    implicit none

    double precision, parameter :: lam = 5.0d0
    double precision, parameter :: angle = 0.25d0 * pi
    double precision, parameter :: evec(2) = (/cos(angle), sin(angle)/)
    double precision :: h, ab, B11, B12, B22, pos(2, 2), error

    n_parcels = 1
    call parcel_alloc(2)

    ab = 1.0d0

    parcels%position(1, :) = 0.0d0
    parcels%volume(1, 1) = ab * pi

    B11 = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2
    B12 = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)
    B22 = lam * sin(angle) ** 2 + 1.0 / lam * cos(angle) ** 2

    parcels%B(1, 1) = B11
    parcels%B(1, 2) = B12

    ! analytic split
    h = 0.25 * sqrt(3.0d0 * ab * lam)
    B11 = 2.0d0 * B11 - 1.5d0 * lam * evec(1) ** 2
    B12 = 2.0d0 * B12 - 1.5d0 * lam * evec(1) * evec(2)
    pos(1, :) = parcels%position(1, :) + h * evec
    pos(2, :) = parcels%position(1, :) - h * evec

    ! numerical split
    call split_ellipses(parcels, threshold=4.0d0)

    !
    ! check result
    !

    error = 0.0d0

    ! first parcel
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, sum(abs(pos(1, :) - parcels%position(1, :))))
    error = max(error, sum(abs(0.5d0 * ab * pi - parcels%volume(1, :))))


    ! second parcel
    error = max(error, abs(parcels%B(2, 1) - B11))
    error = max(error, abs(parcels%B(2, 2) - B12))
    error = max(error, sum(abs(pos(2, :) - parcels%position(2, :))))
    error = max(error, sum(abs(0.5d0 * ab * pi - parcels%volume(2, :))))
    error = max(error, dble(abs(n_parcels - 2)))

    if (error > 1.0e-15) then
        print '(a20, a16)', 'Test ellipse split:', 'FAILED'
    else
        print '(a20, a16)', 'Test ellipse split:', 'PASSED'
    endif

    call parcel_dealloc

end program test_ellipse_split
