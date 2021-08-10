! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks A --> B <--> C <-- D.
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

    call parcel_alloc(4)
    n_parcels = 4

    call parcel_setup((/1, 2, 3, 4/))

    ! geometric merge
    parcel%lambda_max = five
    parcel%vmin_fraction = three

    call find_nearest(isma, ibig, n_merge)

    failed = (n_merge .ne. 1)
    failed = (failed .or. ((isma(1) .ne. 1) .or. (ibig(1) .ne. 2)))

    call print_result_logical('Test nearest algorithm 1', failed)

    contains

        subroutine parcel_setup(ordering)
            integer, intent(in) :: ordering(4)

            parcels%position(ordering(1), 1) = -0.3d0
            parcels%position(ordering(1), 2) = zero
            parcels%volume(ordering(1)) = 0.1d0 * pi

            parcels%position(ordering(2), 1) = -0.1d0
            parcels%position(ordering(2), 2) = zero
            parcels%volume(ordering(2)) = 0.12d0 * pi

            parcels%position(ordering(3), 1) = 0.1d0
            parcels%position(ordering(3), 2) = zero
            parcels%volume(ordering(3)) = 0.12d0 * pi

            parcels%position(ordering(4), 1) = 0.3d0
            parcels%position(ordering(4), 2) = zero
            parcels%volume(ordering(4)) = 0.1d0 * pi
        end subroutine parcel_setup

end program test_nearest_1
