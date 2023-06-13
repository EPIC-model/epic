! =============================================================================
!                       Test ellipsoid split and merge
!
!         This unit test splits an ellispoid and the merges it again. The
!         initial and final ellispoid should be the same. The ellipsoid
!         is centred at the origin and the orientation is varied.
! =============================================================================
program test_ellipsoid_split_merge
    use unit_test
    use constants, only : pi, zero, one, three, four, five, ten, f12, f14, f34
    use parcel_container
    use parcel_ellipsoid, only : get_B33
    use parcel_merge, only : merge_parcels, merge_timer
    use parcel_nearest
    use parcel_split_mod, only : parcel_split, split_timer
    use parameters, only : update_parameters, nx, ny, nz, extent, lower, set_amax, set_vmin, vmin
    use options, only : parcel
    use timer
    implicit none

    double precision, parameter :: lam = five
    double precision :: theta, phi, st, ct, sp, cp
    double precision :: B11, B12, B13, B22, B23, B33, abc, ab
    double precision :: error, a2, b2, c2
    integer :: i, j

    nx = 1
    ny = 1
    nz = 1
    extent = (/ten, ten, ten/)
    lower = (/-five, -five, -five/)
    call update_parameters

    parcel%lambda_max = four


    abc = one
    ab = f34

    a2 = ab * lam
    b2 = ab / lam
    c2 = (abc / ab) ** 2

    ! set vmin to initial parcel volume
    call set_vmin(four / three * abc * pi)

    ! set vmax to half the initial volume
    call set_amax(f12 * a2)

    call register_timer('parcel split', split_timer)
    call register_timer('parcel merge', merge_timer)
    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    call parcel_alloc(2)


    error = zero

    do i = 0, 7
        theta = dble(i) * f14 * pi
        st = dsin(theta)
        ct = dcos(theta)
        do j = 0, 7
            phi = dble(j) * f14 * pi
            sp = dsin(phi)
            cp = dcos(phi)

            call setup_parcels

            ! numerical split
            call parcel_split

            ! make sure we split
            error = max(error, abs(dble(2 - n_parcels)))

            call merge_parcels(parcels)

            ! we should get original parcel
            call check_result
        enddo
    enddo

    call print_result_dp('Test ellipsoid split and split', error, atol=1.0e-14)

    call parcel_dealloc


    contains

        subroutine setup_parcels
            n_parcels = 1

            parcels%position(:, 1) = zero
            parcels%volume(1) = four / three * abc * pi
            parcels%buoyancy(1) = one
#ifndef ENABLE_DRY_MODE
            parcels%humidity(1) = one
#endif
            ! 7 Nov 2021
            ! https://mathworld.wolfram.com/SphericalCoordinates.html
            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp
            B33 = a2 * cp ** 2 + c2 * sp ** 2

            parcels%B(1, 1) = B11
            parcels%B(2, 1) = B12
            parcels%B(3, 1) = B13
            parcels%B(4, 1) = B22
            parcels%B(5, 1) = B23
        end subroutine setup_parcels

        subroutine check_result
            error = max(error, abs(dble(1 - n_parcels)))
            error = max(error, abs(parcels%B(1, 1) - B11))
            error = max(error, abs(parcels%B(2, 1) - B12))
            error = max(error, abs(parcels%B(3, 1) - B13))
            error = max(error, abs(parcels%B(4, 1) - B22))
            error = max(error, abs(parcels%B(5, 1) - B23))
            error = max(error, abs(get_B33(parcels%B(:, 1), parcels%volume(1)) - B33))
            error = max(error, sum(abs(parcels%position(:, 1))))
            error = max(error, abs(four / three * abc * pi - parcels%volume(1)))
            error = max(error, abs(parcels%buoyancy(1) - one))
#ifndef ENABLE_DRY_MODE
            error = max(error, abs(parcels%humidity(1) - one))
#endif
        end subroutine check_result

end program test_ellipsoid_split_merge
