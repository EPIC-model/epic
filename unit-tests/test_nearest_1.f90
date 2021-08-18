! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks:
!               (2a) A(1) <--> B(1)
!               (2b) A(1)  --> B(2)
! =============================================================================
program test_nearest_1
    use unit_test
    use constants, only : pi, zero, two, three, five
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    implicit none

    logical                            :: passed = .true.
    integer, allocatable, dimension(:) :: isma
    integer, allocatable, dimension(:) :: ibig
    integer                            :: n_merge

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call update_parameters

    call parcel_alloc(2)

    !
    !   (2a)
    !
    call parcel_setup_2a

    call find_nearest(isma, ibig, n_merge)

    passed = (n_merge == 1)
    passed = (passed .and. (isma(1) == 1) .and. (ibig(1) == 2))
    call print_result_logical('Test nearest algorithm (setup 2a)', passed)

    !
    !   (2b)
    !
    passed = .true.

    call parcel_setup_2b

    call find_nearest(isma, ibig, n_merge)

    passed = (n_merge == 1)
    passed = (passed .and. (isma(1) == 1) .and. (ibig(1) == 2))
    call print_result_logical('Test nearest algorithm (setup 2b)', passed)

    contains

        subroutine parcel_setup_2a
            n_parcels = 2
            parcels%position(1, 1) = -0.25d0
            parcels%position(1, 2) = zero
            parcels%volume(1) = 0.1d0 * pi

            parcels%position(2, 1) = 0.25d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = 0.1d0 * pi
        end subroutine parcel_setup_2a

        subroutine parcel_setup_2b
            n_parcels = 2
            parcels%position(1, 1) = -0.25d0
            parcels%position(1, 2) = zero
            parcels%volume(1) = 0.1d0 * pi

            parcels%position(2, 1) = 0.25d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = 0.12d0 * pi
        end subroutine parcel_setup_2b


end program test_nearest_1
