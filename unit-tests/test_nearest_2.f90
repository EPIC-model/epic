! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (3a) A(1) <--> B(1) <-- C(1)
!               (3b) A(1) --> B(1) --> C(2)
!               (3c) A(1) --> B(2) <--> C(2)
!               (3d) A(1) --> B(2) <-- C(1)
!               (3e) A(1) --> B(3) <-- C(2)
!               (3f) A(1) --> B(2) --> C(3)
! =============================================================================
program test_nearest_2
    use unit_test
    use permute, only : permute_generate, permute_dealloc, n_permutes, permutes
    use constants, only : pi, zero, two, three, five
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    implicit none

    logical              :: passed = .true.
    integer, allocatable :: isma(:)
    integer, allocatable :: ibig(:)
    integer              :: n_merge, n
    integer              :: ordering(3)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call update_parameters


    call parcel_alloc(3)
    n_parcels = 3

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3a)
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3a(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 2))

        ! B is in ibig since A --> B and C --> B
        passed = (passed .and. &
                    (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))

        ! A and C are in isma
        if (ordering(1) < ordering(3)) then
            passed = (passed .and. &
                        (isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))
        else
            passed = (passed .and. &
                        (isma(1) == ordering(3)) .and. (isma(2) == ordering(1)))
        endif
    enddo

    call print_result_logical('Test nearest algorithm: A(1) <-> B(1) <-  C(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3b)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3b(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 1))

        ! A --> B
        passed = (passed .and. (isma(1) == ordering(1)))
        passed = (passed .and. (ibig(1) == ordering(2)))
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(1)  -> C(2)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3c)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3c(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 2))

        ! B is in ibig since A --> B and C --> B
        passed = (passed .and. &
                    (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))

        ! A and C are in isma
        if (ordering(1) < ordering(3)) then
            passed = (passed .and. &
                        (isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))
        else
            passed = (passed .and. &
                        (isma(1) == ordering(3)) .and. (isma(2) == ordering(1)))
        endif

    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2) <-> C(2)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3d)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3d(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 2))

        ! B is in ibig since A --> B and C --> B
        passed = (passed .and. &
                    (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))

        ! A and C are in isma
        if (ordering(1) < ordering(3)) then
            passed = (passed .and. &
                        (isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))
        else
            passed = (passed .and. &
                        (isma(1) == ordering(3)) .and. (isma(2) == ordering(1)))
        endif

    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2) <-  C(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3e)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3e(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 2))

        ! B is in ibig since A --> B and C --> B
        passed = (passed .and. &
                    (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))

        ! A and C are in isma
        if (ordering(1) < ordering(3)) then
            passed = (passed .and. &
                        (isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))
        else
            passed = (passed .and. &
                        (isma(1) == ordering(3)) .and. (isma(2) == ordering(1)))
        endif

    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(3) <-  C(2)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3f)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3f(ordering)

        call find_nearest(isma, ibig, n_merge)

        passed = (passed .and. (n_merge == 1))

        ! B is in ibig
        passed = (passed .and. (ibig(1) == ordering(2)))

        ! A is in isma
        passed = (passed .and. (isma(1) == ordering(1)))
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2)  -> C(3)', passed)


    call permute_dealloc

    contains

        subroutine parcel_setup_3a(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.1d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.12d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi
        end subroutine parcel_setup_3a

        subroutine parcel_setup_3b(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi
        end subroutine parcel_setup_3b

        subroutine parcel_setup_3c(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi
        end subroutine parcel_setup_3c

        subroutine parcel_setup_3d(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.2d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi
        end subroutine parcel_setup_3d

        subroutine parcel_setup_3e(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.14d0 * pi

            parcels%position(p(3), 1) = 0.2d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi
        end subroutine parcel_setup_3e

        subroutine parcel_setup_3f(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.2d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.14d0 * pi
        end subroutine parcel_setup_3f

end program test_nearest_2
