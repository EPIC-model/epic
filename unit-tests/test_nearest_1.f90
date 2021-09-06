! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (2a) a = b
!               (2b) a - B
! =============================================================================
program test_nearest_1
    use unit_test
    use permute, only : permute_generate, permute_dealloc, n_permutes, permutes
    use constants, only : pi, zero, two, five, ten
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    use timer
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma
    integer, allocatable, dimension(:) :: iclo
    integer                            :: n_merge, n, ordering(2)

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call register_timer('merge nearest', merge_nearest_timer)
    call register_timer('merge tree resolve', merge_tree_resolve_timer)

    parcel%lambda_max = five
    parcel%min_vratio = ten

    call update_parameters

    call parcel_alloc(2)
    n_parcels = 2

    call permute_generate(n_parcels)

    !
    !   (2a) a = b
    !
    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_2a(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AtoB_or_BtoA
    enddo

    call print_result_logical('Test nearest algorithm: a = b', passed)

    !
    !   (2b) a - B
    !
    passed = .true.

    do n = 1, n_permutes
        ordering = permutes(n, :)
        call parcel_setup_2b(ordering)

        call find_nearest(isma, iclo, n_merge)

        call AtoB
    enddo

    call print_result_logical('Test nearest algorithm: a - B', passed)

    contains

        ! (2a) a = b
        subroutine parcel_setup_2a(p)
            integer, intent(in) :: p(2)

            parcels%position(p(1), 1) = -0.25d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.25d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.2d0 * pi
        end subroutine parcel_setup_2a

        ! (2b) a - B
        subroutine parcel_setup_2b(p)
            integer, intent(in) :: p(2)

            parcels%position(p(1), 1) = -0.25d0
            parcels%position(p(1), 2) = zero
            parcels%volume(p(1)) = 0.1d0 * pi

            parcels%position(p(2), 1) = 0.25d0
            parcels%position(p(2), 2) = zero
            parcels%volume(p(2)) = 0.4d0 * pi
        end subroutine parcel_setup_2b

        subroutine AtoB_or_BtoA
            logical :: a_to_b, b_to_a

            passed = (passed .and. (n_merge == 1))

            a_to_b = ((isma(1) == ordering(1)) .and. (iclo(1) == ordering(2)))
            b_to_a = ((iclo(1) == ordering(1)) .and. (isma(1) == ordering(2)))
            passed = (passed .and. (a_to_b .or. b_to_a))
        end subroutine AtoB_or_BtoA

        subroutine AtoB
            logical :: a_to_b
            passed = (passed .and. (n_merge == 1))

            a_to_b = ((isma(1) == ordering(1)) .and. (iclo(1) == ordering(2)))
            passed = (passed .and. a_to_b)
        end subroutine AtoB


end program test_nearest_1
