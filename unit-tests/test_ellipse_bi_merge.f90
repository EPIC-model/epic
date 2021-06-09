! =============================================================================
!                       Test ellipse bi merge
!
!         This unit test checks the merging of two ellipses. The two
!         ellipses, i.e., circles, are located at the origin. The final
!         ellipse is again a circle located at the origin.
! =============================================================================
program test_ellipse_bi_merge
    use unit_test
    use constants, only : pi, one, two, five
    use parcel_container
    use parcel_merge, only : merge_ellipses
    use options, only : parcel, box
    use parameters, only : update_parameters
    use parcel_ellipse
    implicit none

    double precision :: ab, a1b1, a2b2, B11, B12, B22, error, vol, hum, buoy

    box%nc = (/1, 1/)
    box%extent = (/pi, pi/)

    call update_parameters()

    n_parcels = 2
    call parcel_alloc(2)


    ! parcels
    a1b1 = two
    parcels%position(1, :) = zero
    parcels%volume(1) = a1b1 * pi
    parcels%B(1, 1) = a1b1
    parcels%B(1, 2) = zero
    parcels%buoyancy(1) = 1.5d0
    parcels%humidity(1) = 1.3d0

    a2b2 = 0.5d0
    parcels%position(2, :) = zero
    parcels%volume(2) = a2b2 * pi
    parcels%B(2, 1) = a2b2
    parcels%B(2, 2) = zero
    parcels%buoyancy(2) = 1.8d0
    parcels%humidity(2) = 1.2d0

    ! geometric merge
    parcel%lambda = five
    parcel%merge_type = 'bi-geometric'
    parcel%vfraction = two

    call merge_ellipses(parcels)

    ! reference solution
    ab = a1b1 + a2b2  ! a == b since it is a circle
    B11 = ab
    B12 = zero
    B22 = ab
    vol = ab * pi
    buoy = (1.5d0 * a1b1 + 1.8d0 * a2b2) / ab
    hum  = (1.3d0 * a1b1 + 1.2d0 * a2b2) / ab

    !
    ! check result
    !
    error = zero

    error = max(error, abs(dble(n_parcels - 1)))
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, abs(get_B22(parcels%B(1, 1), &
                                   parcels%B(1, 2), &
                                   parcels%volume(1)) - B22))
    error = max(error, sum(abs(parcels%position(1, :))))
    error = max(error, abs(parcels%volume(1) - vol))
    error = max(error, abs(parcels%buoyancy(1) - buoy))
    error = max(error, abs(parcels%humidity(1) - hum))

    call print_result_dp('Test ellipse bi-merge (geometric)', error)

    ! parcels
    n_parcels = 2
    a1b1 = two
    parcels%position(1, :) = zero
    parcels%volume(1) = a1b1 * pi
    parcels%B(1, 1) = a1b1
    parcels%B(1, 2) = zero
    parcels%buoyancy(1) = 1.5d0
    parcels%humidity(1) = 1.3d0

    a2b2 = 0.5d0
    parcels%position(2, :) = zero
    parcels%volume(2) = a2b2 * pi
    parcels%B(2, 1) = a2b2
    parcels%B(2, 2) = zero
    parcels%buoyancy(2) = 1.8d0
    parcels%humidity(2) = 1.2d0

    ! optimal merge
    parcel%merge_type = 'bi-optimal'

    call merge_ellipses(parcels)

    !
    ! check result
    !
    error = zero

    error = max(error, abs(dble(n_parcels - 1)))
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, abs(get_B22(parcels%B(1, 1), &
                                   parcels%B(1, 2), &
                                   parcels%volume(1)) - B22))
    error = max(error, sum(abs(parcels%position(1, :))))
    error = max(error, abs(parcels%volume(1) - vol))
    error = max(error, abs(parcels%buoyancy(1) - buoy))
    error = max(error, abs(parcels%humidity(1) - hum))


    call print_result_dp('Test ellipse bi-merge (optimal)', error)

    call parcel_dealloc

end program test_ellipse_bi_merge
