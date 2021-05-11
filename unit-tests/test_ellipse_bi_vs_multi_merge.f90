! =============================================================================
!                     Test ellipse bi-merge vs multi-merge
!
!         This unit test is like test_ellipse_multi_merge_2 but the ellipses
!         are centred at (1.5, 0.2). Hence, this checks periodicity in x.
! =============================================================================
program test_ellipse_bi_vs_multi_merge
    use unit_test
    use constants, only : pi, one, two, four
    use parcel_container
    use parcel_merge, only : merge_ellipses
    use options, only : parcel_info, box
    use parameters, only : update_parameters
    use ellipse
    implicit none

    double precision :: error, B11, B12, vol, pos(2)

    box%nc = (/1, 1/)
    box%extent = (/pi, pi/)

    call update_parameters()

    call parcel_alloc(2)

    !
    ! muti-geometric merging
    !

    call parcel_setup

    ! geometric merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'multi-geometric'
    parcel_info%vfraction = 3

    call merge_ellipses(parcels)

    B11 = parcels%B(1, 1)
    B12 = parcels%B(1, 2)
    pos = parcels%position(1, :)
    vol = parcels%volume(1, 1)

    call parcel_setup
    parcel_info%merge_type = 'bi-geometric'
    call merge_ellipses(parcels)

    ! check result
    error = 0.0d0
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, sum(abs(parcels%position(1, :) - pos)))
    error = max(error, abs(parcels%volume(1, 1) - vol))

    call print_result_dp('Test ellipse bi-merge vs multi-merge (geometric)', error)

    !
    ! muti-optimal merging
    !

    call parcel_setup

    ! optimal merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'multi-optimal'
    parcel_info%vfraction = 3

    call merge_ellipses(parcels)

    B11 = parcels%B(1, 1)
    B12 = parcels%B(1, 2)
    pos = parcels%position(1, :)
    vol = parcels%volume(1, 1)

    call parcel_setup
    parcel_info%merge_type = 'bi-optimal'
    call merge_ellipses(parcels)

    ! check result
    error = 0.0d0
    error = max(error, abs(parcels%B(1, 1) - B11))
    error = max(error, abs(parcels%B(1, 2) - B12))
    error = max(error, sum(abs(parcels%position(1, :) - pos)))
    error = max(error, abs(parcels%volume(1, 1) - vol))

    call print_result_dp('Test ellipse bi-merge vs multi-merge (optimal)', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            double precision :: a1b1, a2b2, d

            a1b1 = 1.44d0
            a2b2 = 0.25d0

            d = (sqrt(a1b1) + sqrt(a2b2)) * 0.5d0 * sqrt(2.0d0)

            n_parcels = 2
            parcels%position(1, 1) = 1.5d0
            parcels%position(1, 2) = 0.2d0
            parcels%volume(1, 1) = a1b1 * pi
            parcels%B(1, 1) = a1b1
            parcels%B(1, 2) = 0.0d0

            ! small parcel left
            parcels%position(2, 1) = 1.5d0 - d
            parcels%position(2, 2) = 0.2d0 - d
            parcels%volume(2, 1) = a2b2 * pi
            parcels%B(2, 1) = a2b2
            parcels%B(2, 2) = 0.0d0

        end subroutine parcel_setup

end program test_ellipse_bi_vs_multi_merge
