! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (4a) A(1) --> B(1) <--> C(1) <--  D(1)
!               (4b) A(1) --> B(1)  --> C(1)  --> D(2)
!               (4c) A(1) --> B(1)  --> C(2) <--> D(2)
!               (4d) A(1) --> B(2) <--> C(2) <--  D(1)
!               (4e) A(1) --> B(2) <--  C(1) <--  D(1)
!               (4f) A(1) --> B(3) <--  C(2) <--  D(1)
!               (4g) A(1) --> B(2)  --> C(3) <--> D(3)
!               (4h) A(1) --> B(3) <--  C(2) <--  D(1)
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
    integer              :: ordering(4)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call update_parameters


    call parcel_alloc(4)
    n_parcels = 4

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4a) A(1) --> B(1) <--> C(1) <--  D(1)
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4a(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_DtoC
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(1) <-> C(1) <-  D(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4b) A(1) --> B(1)  --> C(1)  --> D(2)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4b(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_CtoD
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(1)  -> C(1)  -> D(2)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4c) A(1) --> B(1)  --> C(2) <--> D(2)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4c(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AB_and_DC
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(1)  -> C(2) <-> D(2)', passed)


    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4d) A(1) --> B(2) <--> C(2) <--  D(1)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4d(ordering)

        call find_nearest(isma, ibig, n_merge)

!         print *, "--------------------"
!         print *, ordering
!         print *, n_merge
!         print *, isma(1:n_merge)
!         print *, ibig(1:n_merge)
        call AtoB_and_DtoC
    enddo
!     stop

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2) <-> C(2) <-  D(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4e) A(1) --> B(2) <--  C(1) <--  D(1)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4e(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_DtoC
!         call ACtoB
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2) <-  C(1) <-  D(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4f) A(1) --> B(3) <--  C(2) <--  D(1)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4f(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_DtoC
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(3) <-  C(2) <-  D(1)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4g) A(1) --> B(2)  --> C(3) <--> D(3)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4g(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_DtoC
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2)  -> C(3) <-> D(3)', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (4h) A(1) --> B(3) <--  C(2) <--  D(1)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_4h(ordering)

        call find_nearest(isma, ibig, n_merge)

        call AtoB_and_DtoC
    enddo

    call print_result_logical('Test nearest algorithm: A(1)  -> B(2)  -> C(3) <-> D(3)', passed)

    call permute_dealloc

    contains

        ! (4a) A(1) --> B(1) <--> C(1) <--  D(1)
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

        ! (4b) A(1) --> B(1)  --> C(1)  --> D(2)
        subroutine parcel_setup_4b(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.1d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.15d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.12d0 * pi
        end subroutine parcel_setup_4b

        ! (4c) A(1) --> B(1)  --> C(2) <--> D(2)
        subroutine parcel_setup_4c(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.25d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.12d0 * pi
        end subroutine parcel_setup_4c

        ! (4d) A(1) --> B(2) <--> C(2) <--  D(1)
        subroutine parcel_setup_4d(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.35d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.1d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi

            parcels%position(p(4), 1) = 0.35d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4d

        ! (4e) A(1) --> B(2) <--  C(1) <--  D(1)
        subroutine parcel_setup_4e(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.1d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.0d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.15d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4e

        ! (4f) A(1) --> B(3) <--  C(2) <--  D(1)
        subroutine parcel_setup_4f(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.2d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.1d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.14d0 * pi

            parcels%position(p(3), 1) = 0.0d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi

            parcels%position(p(4), 1) = 0.1d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4f

        ! (4g) A(1) --> B(2)  --> C(3) <--> D(3)
        subroutine parcel_setup_4g(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.15d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.12d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.14d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.14d0 * pi
        end subroutine parcel_setup_4g

        ! (4h) A(1) --> B(3) <--  C(2) <--  D(1)
        subroutine parcel_setup_4h(p)
            integer, intent(in) :: p(4)
            parcels%position(p(1), 1) = -0.15d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.05d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.14d0 * pi

            parcels%position(p(3), 1) = 0.1d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.12d0 * pi

            parcels%position(p(4), 1) = 0.2d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi
        end subroutine parcel_setup_4h

!         subroutine BtoC_or_CtoB
!             passed = (passed .and. (n_merge == 1))
!
!             if (ordering(2) < ordering(3)) then
!                 ! B --> C
!                 passed = (passed .and. (isma(1) == ordering(2)))
!                 passed = (passed .and. (ibig(1) == ordering(3)))
!             else
!                 ! C --> B
!                 passed = (passed .and. (isma(1) == ordering(3)))
!                 passed = (passed .and. (ibig(1) == ordering(2)))
!             endif
!         end subroutine BtoC_or_CtoB

        subroutine AB_and_DC
            passed = (passed .and. (n_merge == 2))

            ! A and B
            passed = (passed .and.                                                       &
                        (((isma(1) == ordering(1)) .and. (ibig(1) == ordering(2)))  .or. &
                         ((ibig(1) == ordering(1)) .and. (isma(1) == ordering(2))))      &
                    )
            ! C and D
            passed = (passed .and.                                                       &
                        (((isma(2) == ordering(3)) .and. (ibig(2) == ordering(4)))  .or. &
                         ((ibig(2) == ordering(3)) .and. (isma(2) == ordering(4))))      &
                    )

            print *, passed
        end subroutine AB_and_DC

        subroutine AtoB_and_DtoC
            passed = (passed .and. (n_merge == 2))

            ! A and D are in isma
            passed = (passed .and.                                                       &
                        (((isma(1) == ordering(1)) .and. (isma(2) == ordering(4)))  .or. &
                         ((isma(2) == ordering(1)) .and. (isma(1) == ordering(4))))      &
                    )
            ! B and C are in ibig
            passed = (passed .and.                                                       &
                        (((ibig(1) == ordering(2)) .and. (ibig(2) == ordering(3)))  .or. &
                         ((ibig(2) == ordering(2)) .and. (ibig(1) == ordering(3))))      &
                    )
        end subroutine AtoB_and_DtoC

        subroutine AtoB_and_CtoD
            passed = (passed .and. (n_merge == 2))

            ! A and C are in isma
            passed = (passed .and.                                                       &
                        (((isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))  .or. &
                         ((isma(2) == ordering(1)) .and. (isma(1) == ordering(3))))      &
                    )
            ! B and D are in ibig
            passed = (passed .and.                                                       &
                        (((ibig(1) == ordering(2)) .and. (ibig(2) == ordering(4)))  .or. &
                         ((ibig(2) == ordering(2)) .and. (ibig(1) == ordering(4))))      &
                    )
        end subroutine AtoB_and_CtoD

!         subroutine BDtoC
!             passed = (passed .and. (n_merge == 2))
!
!             ! C is in ibig since B --> C and D --> C
!             passed = (passed .and. &
!                         (ibig(1) == ordering(3)) .and. (ibig(2) == ordering(3)))
!
!             ! B and D are in isma
!             if (ordering(2) < ordering(4)) then
!                 passed = (passed .and. &
!                             (isma(1) == ordering(2)) .and. (isma(2) == ordering(4)))
!             else
!                 passed = (passed .and. &
!                             (isma(1) == ordering(4)) .and. (isma(2) == ordering(2)))
!             endif
!         end subroutine BDtoC

!         subroutine ACtoB
!             passed = (passed .and. (n_merge == 2))
!
!             ! B is in ibig since A --> B and C --> B
!             passed = (passed .and. &
!                         (ibig(1) == ordering(2)) .and. (ibig(2) == ordering(2)))
!
!             ! A and C are in isma
!             if (ordering(1) < ordering(3)) then
!                 passed = (passed .and. &
!                             (isma(1) == ordering(1)) .and. (isma(2) == ordering(3)))
!             else
!                 passed = (passed .and. &
!                             (isma(1) == ordering(3)) .and. (isma(2) == ordering(1)))
!             endif
!         end subroutine ACtoB

!         subroutine AtoB
!             passed = (passed .and. (n_merge == 1))
!
!             ! A is in isma
!             passed = (passed .and. (isma(1) == ordering(1)))
!
!             ! B is in ibig
!             passed = (passed .and. (ibig(1) == ordering(2)))
!         end subroutine AtoB

end program test_nearest_2
