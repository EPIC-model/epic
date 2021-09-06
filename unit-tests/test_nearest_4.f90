! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (5a) (a, b) - c - d = e
!               (5b) (a, b) - C - (d, e)
! =============================================================================
program test_nearest_4
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
    integer              :: ordering(5)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters

    call parcel_alloc(5)
    n_parcels = 5

    call permute_generate(n_parcels)

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !   (5a) (a, b) - c - d = e
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_5a(ordering)

        call find_nearest(isma, iclo, n_merge)

        call ABtoC_and_DE
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

        call find_nearest(isma, iclo, n_merge)

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

        subroutine ABtoC_and_DE
            integer :: ia,ib,ide,ii

            passed = (passed .and. (n_merge == 3))

            do ii=1,3
              if(isma(ii)==ordering(1)) then
                ia=ii
              elseif(isma(ii)==ordering(2)) then
                ib=ii
              elseif(isma(ii)==ordering(4) .or. (isma(ii)==ordering(5) )) then
                ide=ii
              endif
            enddo

            ! A to C
            passed = (passed .and.                                                         &
                        ((isma(ia) == ordering(1)) .and. (iclo(ia) == ordering(3)))        &
                    )

            ! B to C
            passed = (passed .and.                                                         &
                        ((isma(ib) == ordering(2)) .and. (iclo(ib) == ordering(3)))        &
                    )

            ! D and E
            passed = (passed .and.                                                         &
                        (((isma(ide) == ordering(4)) .and. (iclo(ide) == ordering(5)))  .or. &
                         ((iclo(ide) == ordering(4)) .and. (isma(ide) == ordering(5))))      &
                    )

        end subroutine ABtoC_and_DE

        subroutine ABDEtoC
            integer :: i
            passed = (passed .and. (n_merge == 4))

            ! C is the only big parcel
            do i = 1, 4
                passed = (passed .and. (iclo(i) == ordering(3)))
            enddo
        end subroutine ABDEtoC

end program test_nearest_4
