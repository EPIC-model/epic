! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks A <--> B (equal-sized parcels).
! =============================================================================
program test_nearest_1
    use unit_test
    use constants, only : pi, zero, two, three, five
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    implicit none

    logical                            :: failed = .false.
    integer, allocatable, dimension(:) :: isma
    integer, allocatable, dimension(:) :: ibig
    integer                            :: n_merge

    nx = 1
    nz = 1
    lower  = (/-pi / two, -pi /two/)
    extent = (/pi, pi/)

    call update_parameters

    call parcel_alloc(2)

    call parcel_setup

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call find_nearest(isma, ibig, n_merge)

    failed = (n_merge .ne. 1)
    failed = (failed .or. ((isma(1) .ne. 1) .or. (ibig(1) .ne. 2)))

    call print_result_logical('Test nearest algorithm 1', failed)

    contains

        subroutine parcel_setup
            n_parcels = 2
            parcels%position(1, 1) = -0.25d0
            parcels%position(1, 2) = zero
            parcels%volume(1) = 0.1d0 * pi

            parcels%position(2, 1) = 0.25d0
            parcels%position(2, 2) = zero
            parcels%volume(2) = 0.1d0 * pi
        end subroutine parcel_setup


end program test_nearest_1
