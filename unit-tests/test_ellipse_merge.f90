! =============================================================================
!                       Test ellipse bi merge
!
!         This unit test checks the merging of two ellipses. The two
!         ellipses, i.e., circles, are located at the origin. The final
!         ellipse is again a circle located at the origin.
! =============================================================================
program test_ellipse_bi_merge
    use constants, only : pi
    use parcel_container
    use parcel_merge, only : merge_ellipses
    use options, only : parcel_info, grid
    use parameters, only : update_parameters
    use ellipse
    implicit none

    double precision :: ab, a1b1, a2b2, B11, B12, B22, error, vol

    grid = (/21, 21/)

    call update_parameters()

    n_parcels = 2
    call parcel_alloc(2)


    ! parcels
    a1b1 = 2.0d0
    parcels%position(1, :) = 0.0d0
    parcels%volume(1, 1) = a1b1 * pi
    parcels%B(1, 1) = a1b1
    parcels%B(1, 2) = 0.0d0

    a2b2 = 0.5d0
    parcels%position(2, :) = 0.0d0
    parcels%volume(2, 1) = a2b2 * pi
    parcels%B(2, 1) = a2b2
    parcels%B(2, 2) = 0.0d0

    ! geometric merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'bi-geometric'
    parcel_info%vfraction = 1.0e-2

    call merge_ellipses(parcels)

    ! reference solution
    ab = a1b1 + a2b2  ! a == b since it is a circle
    B11 = ab
    B12 = 0.0d0
    B22 = ab
    vol = ab * pi

    !
    ! check result
    !
    error = 0.0d0

    error = max(error, abs(dble(n_parcels - 1)))
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, abs(get_B22(parcels%B(1, 1), &
                                   parcels%B(1, 2), &
                                   parcels%volume(1, 1)) - B22))
    error = max(error, sum(abs(parcels%position(1, :))))
    error = max(error, abs(parcels%volume(1, 1) - vol))

    if (error > 1.0e-15) then
        print '(a26, a10)', 'Test ellipse merge (geo):', 'FAILED'
    else
        print '(a26, a10)', 'Test ellipse merge (geo):', 'PASSED'
    endif


    ! parcels
    n_parcels = 2
    a1b1 = 2.0d0
    parcels%position(1, :) = 0.0d0
    parcels%volume(1, 1) = a1b1 * pi
    parcels%B(1, 1) = a1b1
    parcels%B(1, 2) = 0.0d0

    a2b2 = 0.5d0
    parcels%position(2, :) = 0.0d0
    parcels%volume(2, 1) = a2b2 * pi
    parcels%B(2, 1) = a2b2
    parcels%B(2, 2) = 0.0d0

    ! optimal merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'bi-optimal'
    parcel_info%vfraction = 1.0e-2

    call merge_ellipses(parcels)

    !
    ! check result
    !
    error = 0.0d0

    error = max(error, abs(dble(n_parcels - 1)))
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, abs(get_B22(parcels%B(1, 1), &
                                   parcels%B(1, 2), &
                                   parcels%volume(1, 1)) - B22))
    error = max(error, sum(abs(parcels%position(1, :))))
    error = max(error, abs(parcels%volume(1, 1) - vol))


    if (error > 1.0e-15) then
        print '(a26, a10)', 'Test ellipse merge (opt):', 'FAILED'
    else
        print '(a26, a10)', 'Test ellipse merge (opt):', 'PASSED'
    endif

    call parcel_dealloc

end program test_ellipse_bi_merge
