! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (4a) a - b = c - d
!               (4b) a - b - c = d
!               (4c) a - B - c - d
!               (4d) a = b   c = d
!               (4e) a - b - c - D
! =============================================================================
program test_nearest_3
    use unit_test
    use permute, only : permute_generate, permute_dealloc, n_permutes, permutes
    use constants, only : pi, zero, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    implicit none

    logical              :: passed = .true.
    integer, allocatable :: isma(:)
    integer, allocatable :: iclo(:)
    integer              :: n_merge, n
    integer              :: ordering(4)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters

    call parcel_alloc(4)
    n_parcels = 4

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4a) a - b = c - d
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4a(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: a - b = c - d', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4b) a - b - c = d
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4b(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: a - b - c = d', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4c) a - B - c - d
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4c(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: a - B - c - d', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4d) a = b   c = d
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4d(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: a = b   c = d', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4e) a - b - c - D
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4e(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: a - b - c - D', passed)

    call permute_dealloc

    contains

        ! (4a) a - b = c - d
        subroutine parcel_setup_4a(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.05d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4a

        ! (4b) a - b - c = d
        subroutine parcel_setup_4b(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.25d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4b

        ! (4c) a - B - c - d
        subroutine parcel_setup_4c(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.1d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.4d0 * pi

            parcels%position(p(3), 1) = 0.0d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi

            parcels%position(p(4), 1) = 0.15d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4c

        ! (4d) a = b   c = d
        subroutine parcel_setup_4d(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.1d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.101d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4d

        ! (4e) a - b - c - D
        subroutine parcel_setup_4e(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.11d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.02d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.05d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.1d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.4d0 * pi
        end subroutine parcel_setup_4e
!
        subroutine AB_and_DC
            integer :: ia, ic   ! order of merge can change

            passed = (passed .and. (n_merge == 2))

            ia = 1
            ic = 2

            if ((isma(2) == ordering(1)) .or. (iclo(2) == ordering(1))) then
                ia = 2
                ic = 1
            endif

            ! A and B
            passed = (passed .and.                                                         &
                        (((isma(ia) == ordering(1)) .and. (iclo(ia) == ordering(2)))  .or. &
                         ((iclo(ia) == ordering(1)) .and. (isma(ia) == ordering(2))))      &
                    )
            ! C and D
            passed = (passed .and.                                                         &
                        (((isma(ic) == ordering(3)) .and. (iclo(ic) == ordering(4)))  .or. &
                         ((iclo(ic) == ordering(3)) .and. (isma(ic) == ordering(4))))      &
                    )

        end subroutine AB_and_DC

end program test_nearest_3
