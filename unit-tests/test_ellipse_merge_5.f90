! =============================================================================
!                       Test ellipse multi merge
!
!       This unit test checks group merging. It places three parcels
!       in a row where the first two parcels are of equal size. In each call
!       to the merging routine
!       only two parcels should merge. After two calls, only one big parcel
!       exists.
! =============================================================================
program test_ellipse_multi_merge_5
    use unit_test
    use constants, only : pi, one, two, three, four, five
    use parcel_container
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_ellipse
    use timer
    implicit none

    double precision :: error

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('parcel merge', merge_timer)

    call update_parameters

    call parcel_alloc(3)

    call parcel_setup

    ! geometric merge
    parcel%lambda_max = five
    parcel%merge_type = 'geometric'
    parcel%vmin_fraction = three

    ! first merge
    call merge_ellipses(parcels)

    ! check result
    error = dble(abs(n_parcels - 2))

    ! check if the small and big parcels merged first
    error = max(error, dabs(parcels%volume(1) - 0.25d0 * pi))

    ! second merge
    call merge_ellipses(parcels)

    ! check result
    error = max(error, dble(abs(n_parcels - 1)))

    call print_result_dp('Test ellipse merge 5', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            n_parcels = 3
            parcels%position(1, 1) = -0.5d0
            parcels%position(1, 2) = zero
            parcels%volume(1) = 0.25d0 * pi
            parcels%B(1, 1) = 0.25d0
            parcels%B(1, 2) = zero

            parcels%position(2, 1) = 0.0d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = 0.25d0 * pi
            parcels%B(2, 1) = 0.25d0
            parcels%B(2, 2) = zero

            parcels%position(3, 1) = 0.5d0
            parcels%position(3, 2) = zero
            parcels%volume(3) = 0.27d0 * pi
            parcels%B(3, 1) = 0.27d0
            parcels%B(3, 2) = zero

        end subroutine parcel_setup

end program test_ellipse_multi_merge_5
