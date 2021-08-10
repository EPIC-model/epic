! =============================================================================
!                       Test nearest algorithm
!
!           This unit test checks A(1) --> B(2) <--> C(2) <-- D(1).
!           [behaves like (2b) and (3c)
! =============================================================================
program test_nearest_4c
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
    integer                            :: permutation(24, 4), i, order(4)

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

    permutation( 1, :) = (/1, 2, 3, 4/)
    permutation( 2, :) = (/1, 2, 4, 3/)
    permutation( 3, :) = (/1, 3, 2, 4/)
    permutation( 4, :) = (/1, 4, 2, 3/)
    permutation( 5, :) = (/1, 3, 4, 2/)
    permutation( 6, :) = (/1, 4, 3, 2/)
    permutation( 7, :) = (/2, 1, 3, 4/)
    permutation( 8, :) = (/2, 1, 4, 3/)
    permutation( 9, :) = (/3, 1, 2, 4/)
    permutation(10, :) = (/3, 1, 4, 2/)
    permutation(11, :) = (/4, 1, 2, 3/)
    permutation(12, :) = (/4, 1, 3, 2/)
    permutation(13, :) = (/2, 3, 1, 4/)
    permutation(14, :) = (/2, 4, 1, 3/)
    permutation(15, :) = (/3, 2, 1, 4/)
    permutation(16, :) = (/4, 2, 1, 3/)
    permutation(17, :) = (/3, 4, 1, 2/)
    permutation(18, :) = (/4, 3, 1, 2/)
    permutation(19, :) = (/2, 3, 4, 1/)
    permutation(20, :) = (/2, 4, 3, 1/)
    permutation(21, :) = (/3, 2, 4, 1/)
    permutation(22, :) = (/4, 2, 3, 1/)
    permutation(23, :) = (/3, 4, 2, 1/)
    permutation(24, :) = (/4, 3, 2, 1/)

    do i = 1, 24
        order = permutation(i, :)
        call parcel_setup(order)
        call find_nearest(isma, ibig, n_merge)

        failed = (failed .or. (n_merge .ne. 2))

        if (isma(1) == order(1)) then
            failed = (failed .or. (ibig(1) .ne. order(2)))
            failed = (failed .or. (ibig(2) .ne. order(3)))
        else
            failed = (failed .or. (ibig(2) .ne. order(2)))
            failed = (failed .or. (ibig(1) .ne. order(3)))
        endif
    enddo

    call print_result_logical('Test nearest algorithm 4c', failed)

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

end program test_nearest_4c
