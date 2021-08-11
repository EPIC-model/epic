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
    use permute
    use constants, only : pi, zero, two, three, five
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, nz
    use parcel_nearest
    implicit none

    logical              :: failed = .false.
    integer, allocatable :: isma(:)
    integer, allocatable :: ibig(:)
    integer              :: n_merge, np, n
    integer, allocatable :: permutes(:, :)
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

    call permute_generate(permutes, n_parcels)

    np = size(permutes) / n_parcels

    !
    !   (3a)
    !
    do n = 1, np
        ordering = permutes(n, :)
        call parcel_setup_3a(ordering)

        call find_nearest(isma, ibig, n_merge)

        failed = (failed .or. (n_merge .ne. 2))

        ! A --> B
        failed = (failed .or. ((isma(1) .ne. ordering(1)) .or. (ibig(1) .ne. ordering(2))))

        ! C --> B
        failed = (failed .or. ((isma(2) .ne. ordering(3)) .or. (ibig(2) .ne. ordering(2))))

        print *, ordering, isma(1), ordering(1)
!         stop
    enddo

    call print_result_logical('Test nearest algorithm: A(1) <-> B(1) <- C(1)', failed)

    stop

    failed = .false.

    !
    !   (3b)
    !
!     call parcel_setup_3b

    call find_nearest(isma, ibig, n_merge)

    failed = (n_merge .ne. 1)
    failed = (failed .or. ((isma(1) .ne. 1) .or. (ibig(1) .ne. 2)))
    call print_result_logical('Test nearest algorithm: A(1)  -> B(1) -> C(2)', failed)

    deallocate(permutes)

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
!
!         subroutine parcel_setup_3b
!             integer, intent(in) :: p(3)
!             parcels%position(p(1), 1) = -0.2d0
!             parcels%position(p(1), 2) = zero
!             parcels%volume(p(1)) = 0.1d0 * pi
!
!             parcels%position(p(2), 1) = 0.0d0
!             parcels%position(p(2), 2) = zero
!             parcels%volume(p(2)) = 0.1d0 * pi
!
!             parcels%position(p(3), 1) = 0.1d0
!             parcels%position(p(3), 2) = zero
!             parcels%volume(p(3)) = 0.12d0 * pi
!         end subroutine parcel_setup_3b


end program test_nearest_2
