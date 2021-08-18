! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (5a) (a, b) - c - d = e
!               (5b) (a, b) - C - (d, e)
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
    integer, allocatable :: ibig(:)
    integer              :: n_merge, n
    integer              :: ordering(5)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call update_parameters


    call parcel_alloc(5)
    n_parcels = 5

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = ten

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (5a) (a, b) - c - d = e
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_5a(ordering)

        call find_nearest(isma, ibig, n_merge)

        call DtoE_or_EtoD
    enddo

    call print_result_logical('Test nearest algorithm: (a, b) - c - d = e', passed)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (5b) (a, b) - C - (d, e)
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_5b(ordering)

        call find_nearest(isma, ibig, n_merge)

        call ABDEtoC
    enddo

    call print_result_logical('Test nearest algorithm: (a, b) - C - (d, e)', passed)

    call permute_dealloc

    contains

        ! (5a) (a, b) - c - d = e
        subroutine parcel_setup_5a(p)
            integer, intent(in) :: p(5)
            parcels%position(p(1), 1) = 0.12d0
            parcels%position(p(1), 2) = 0.07d0
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.12d0
            parcels%position(p(2), 2) = -0.07d0
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.19d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.1d0 * pi

            parcels%position(p(4), 1) = 0.25d0
            parcels%position(p(4), 2) = zero
            parcels%volume(p(4)) = 0.1d0 * pi

            parcels%position(p(5), 1) = 0.3d0
            parcels%position(p(5), 2) = zero
            parcels%volume(p(5)) = 0.1d0 * pi
        end subroutine parcel_setup_5a

        ! (5b) (a, b) - C - (d, e)
        subroutine parcel_setup_5b(p)
            integer, intent(in) :: p(5)
            parcels%position(p(1), 1) = -0.12d0
            parcels%position(p(1), 2) = 0.07d0
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = -0.12d0
            parcels%position(p(2), 2) = -0.07d0
            parcels%volume(p(2)) = 0.1d0 * pi

            parcels%position(p(3), 1) = 0.0d0
            parcels%position(p(3), 2) = zero
            parcels%volume(p(3)) = 0.4d0 * pi

            parcels%position(p(4), 1) = 0.12d0
            parcels%position(p(4), 2) = 0.07d0
            parcels%volume(p(4)) = 0.1d0 * pi

            parcels%position(p(5), 1) = 0.12d0
            parcels%position(p(5), 2) = -0.07d0
            parcels%volume(p(5)) = 0.1d0 * pi
        end subroutine parcel_setup_5b

        subroutine DtoE_or_EtoD
            passed = (passed .and. (n_merge == 1))

            if (ordering(4) < ordering(5)) then
                ! D --> E
                passed = (passed .and. (isma(1) == ordering(4)))
                passed = (passed .and. (ibig(1) == ordering(5)))
            else
                ! E --> D
                passed = (passed .and. (isma(1) == ordering(5)))
                passed = (passed .and. (ibig(1) == ordering(4)))
            endif
        end subroutine DtoE_or_EtoD

        subroutine ABDEtoC
            integer :: i
            passed = (passed .and. (n_merge == 4))

            ! C is the only big parcel
            do i = 1, 4
                passed = (passed .and. (ibig(i) == ordering(3)))
            enddo
        end subroutine ABDEtoC

end program test_nearest_3
