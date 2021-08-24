! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the symmetry by performing a mirrored
!         group-merge (mirror axis x).
! =============================================================================
program test_ellipse_multi_merge_symmetry
    use unit_test
    use constants, only : pi, one, two, three, four, five
    use parcel_container
    use parcel_interpl, only : vol2grid_symmetry_error, sym_vol2grid_timer
    use parcel_merge, only : merge_ellipses, merge_timer
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use fields, only : sym_volg
    use parcel_ellipse
    use timer
    implicit none

    double precision :: error

    nx = 2
    nz = 2
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('parcel merge', merge_timer)
    call register_timer('symmetric vol2grid', sym_vol2grid_timer)

    parcel%lambda_max = five
    parcel%min_vratio = three

    call update_parameters

    allocate(sym_volg(-1:nz+1, 0:nx-1))

    call parcel_alloc(6)

    call parcel_setup

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse group-merge symmetry', error)

    call parcel_dealloc

    deallocate(sym_volg)

    contains

        subroutine parcel_setup
            double precision :: a1b1, a2b2
            integer :: n

            a1b1 = 1.44d0
            a2b2 = 0.25d0

            n_parcels = 6
            parcels%position(1, 1) = -f12
            parcels%position(1, 2) = 0.2d0
            parcels%volume(1) = a1b1 * pi
            parcels%B(1, 1) = 1.2d0 * a1b1
            parcels%B(1, 2) = -0.4d0

            parcels%position(2, 1) = -0.6d0
            parcels%position(2, 2) = 0.3d0
            parcels%volume(2) = a2b2 * pi
            parcels%B(2, 1) = 0.8d0 * a2b2
            parcels%B(2, 2) = f12

            parcels%position(3, 1) = -0.3d0
            parcels%position(3, 2) = -0.1d0
            parcels%volume(3) = a2b2 * pi
            parcels%B(3, 1) = 0.9d0 * a2b2
            parcels%B(3, 2) = -0.1d0

            !
            ! mirrored parcels
            !
            do n = 1, 3
                parcels%position(3+n, 1) = -parcels%position(n, 1)
                parcels%position(3+n, 2) =  parcels%position(n, 2)
                parcels%volume(3+n) = parcels%volume(n)
                parcels%B(3+n, 1) =  parcels%B(n, 1)
                parcels%B(3+n, 2) = -parcels%B(n, 2)
            enddo

        end subroutine parcel_setup

        function eval_max_error() result(max_err)
            double precision :: max_err

            max_err = zero
            max_err = max(max_err, abs(dble(n_parcels - 2)))

            call vol2grid_symmetry_error

            max_err = max(max_err, maxval(abs(sym_volg)))

        end function eval_max_error

end program test_ellipse_multi_merge_symmetry
