! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (3a) a = b - c
!               (3b) a - b - C
!               (3c) a = b   C
!               (3d) a - B - c
! =============================================================================
program test_nearest_2
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
    integer, allocatable :: ibig(:)
    integer              :: n_merge, n
    integer              :: ordering(3)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters


    call parcel_alloc(3)
    n_parcels = 3

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3a) a = b - c
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3a(ordering)

        call find_nearest(isma, ibig, n_merge)

        call ACtoB
    enddo

    call print_result_logical('Test nearest algorithm: a = b - c', passed)


    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3b) a - b - C
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3b(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB
    enddo

    call print_result_logical('Test nearest algorithm: a - b - C', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3c) a = b   C
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3c(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_or_BtoA
    enddo

    call print_result_logical('Test nearest algorithm: a = b   C', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (3d) a - B - c
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_3d(ordering)

        call find_nearest(isma, ibig, n_merge)

        call ACtoB
    enddo

    call print_result_logical('Test nearest algorithm: a - B - c', passed)

    call permute_dealloc

    contains

        ! (3a) a = b - c
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

        ! (3b) a - b - C
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
            parcels%volume(p(3)) = 0.4d0 * pi
        end subroutine parcel_setup_3b

        ! (3c) a = b   C
        subroutine parcel_setup_3c(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.1d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.11d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.4d0 * pi
        end subroutine parcel_setup_3c

        ! (3d) a - B - c
        subroutine parcel_setup_3d(p)
            integer, intent(in) :: p(3)
            parcels%position(p(1), 1) = -0.1d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.0d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.4d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi
        end subroutine parcel_setup_3d

        subroutine AtoB
            passed = (passed .and. (n_merge == 1))

           ! A --> B
           passed = (passed .and. (isma(1) == ordering(1)))
           passed = (passed .and. (ibig(1) == ordering(2)))

        end subroutine AtoB

        subroutine AtoB_or_BtoA
            logical :: a_to_b, b_to_a

            passed = (passed .and. (n_merge == 1))

            a_to_b = ((isma(1) == ordering(1)) .and. (ibig(1) == ordering(2)))
            b_to_a = ((ibig(1) == ordering(1)) .and. (isma(1) == ordering(2)))
            passed = (passed .and. (a_to_b .or. b_to_a))
        end subroutine AtoB_or_BtoA

        subroutine BtoC
            passed = (passed .and. (n_merge == 1))

            ! B --> C
            passed = (passed .and. (isma(1) == ordering(2)))
            passed = (passed .and. (ibig(1) == ordering(3)))
        end subroutine BtoC

        subroutine ACtoB
            passed = (passed .and. (n_merge == 2))

            ! B is in ibig since A --> B and C --> B
            passed = (passed .and. &
                        (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))


            ! A and C are in isma
            passed = (passed .and. &
                        (((isma(1) == ordering(1)) .and. (isma(2) == ordering(3))) .or. &
                         ((isma(1) == ordering(3)) .and. (isma(2) == ordering(1))))     &
                    )
        end subroutine ACtoB

end program test_nearest_2
