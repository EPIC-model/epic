! =============================================================================
!                         Test MPI parcel split
!
!   This unit test checks parcel splitting. Each MPI rank places parcels
!   in each interior cell attached to a halo cell. In each split, one child
!   stays in the domain of *this* process and the other chlid leaves the
!   domain of *this* process.
! =============================================================================
program test_mpi_parcel_split
    use constants, only : zero, one, f18, f14, f12, f34, fpi, pi, four
    use unit_test
    use mpi_environment
    use mpi_layout
    use parcels_mod, only : parcels
    use parcel_container, only : resize_timer
    use parcel_mpi
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, dx, vcell, set_amax
    use mpi_collectives
    use parcel_split_mod, only : parcel_split, split_timer
    use options, only : parcel
    use mpi_timer
    implicit none

    logical          :: passed = .true.
    double precision :: shift = 0.0001d0
    double precision :: abc, a2, b2, c2
    integer          :: i, j, n, n_total
    integer          :: n_orig_total
    integer          :: n_orig_local

    call mpi_env_initialise

    call register_timer('parcel split', split_timer)
    call register_timer('parcel container resize', resize_timer)

    passed = (world%err == 0)

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent = four

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters


    parcels%local_num = 2 * (box%hi(2) - box%lo(2) + 1) + 2 * (box%hi(1) - box%lo(1) + 1)
    n_total = 2 * parcels%local_num

    call parcels%allocate(n_total)

    parcel%lambda_max = four

    n = 1

    ! volume per parcel is f12 * vcell
    ! f12 * vcell = four / three * abc * pi --> abc = f34 * f12 * vcell * fpi
    abc = f34 * f12 * vcell * fpi
    a2 = f34 * abc
    b2 = f18 * abc
    c2 = b2

    abc = sqrt(a2)
    call set_amax(abc)

    ! place parcels in the last interior cells in the west
    i = box%lo(1)
    do j = box%lo(2), box%hi(2)
        parcels%position(1, n) = dble(i) * dx(1) + shift
        parcels%position(2, n) = (dble(j) + f12) * dx(2)
        parcels%position(3, n) = f12
        call fill_shape(n, angle=zero)
        n = n + 1
    enddo

    ! place parcels in the last interior cells in the east
    i = box%hi(1) + 1
    do j = box%lo(2), box%hi(2)
        parcels%position(1, n) = dble(i) * dx(1) - shift
        parcels%position(2, n) = (dble(j) + f12) * dx(2)
        parcels%position(3, n) = f12
        call fill_shape(n, angle=zero)
        n = n + 1
    enddo


    ! place parcels in the last interior cells in the south
    j = box%lo(2)
    do i = box%lo(1), box%hi(1)
        parcels%position(1, n) = (dble(i) + f12) * dx(1)
        parcels%position(2, n) = dble(j) * dx(2) + shift
        parcels%position(3, n) = f12
        call fill_shape(n, angle=f12 * pi)
        n = n + 1
    enddo

    ! place parcels in the last interior cells in the north
    j = box%hi(2) + 1
    do i = box%lo(1), box%hi(1)
        parcels%position(1, n) = (dble(i) + f12) * dx(1)
        parcels%position(2, n) = dble(j) * dx(2) - shift
        parcels%position(3, n) = f12
        call fill_shape(n, angle=f12 * pi)
        n = n + 1
    enddo

    n_orig_total = parcels%local_num
    n_orig_local = parcels%local_num

    call mpi_blocking_reduce(n_orig_total, MPI_SUM, world)

    parcels%volume(1:parcels%local_num) = f12 * vcell
    parcels%vorticity(:, 1:parcels%local_num) = world%rank + 1
    parcels%buoyancy(1:parcels%local_num) = world%rank + 1

    call parcel_split

    passed = (passed .and. (2 * n_orig_local == parcels%local_num))

    if (world%rank == world%root) then
        passed = (passed .and. (2 * n_orig_total == parcels%total_num))
    endif

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    passed = (passed .and. (world%err == 0))

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel split', passed)
    endif

    call mpi_env_finalise

    contains

        subroutine fill_shape(m, angle)
            integer,          intent(in) :: m
            double precision, intent(in) :: angle
            double precision             :: B11, B12, B13, B22, B23, st, ct, sp, cp

            st = sin(angle)
            ct = cos(angle)
            sp = zero ! = sin(phi) with phi = 0
            cp = one  ! = cos(phi) with phi = 0

            B11 = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
            B12 = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
            B13 = (a2 - c2) * ct * sp * cp
            B22 = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
            B23 = (a2 - c2) * st * sp * cp

            parcels%B(1, m) = B11
            parcels%B(2, m) = B12
            parcels%B(3, m) = B13
            parcels%B(4, m) = B22
            parcels%B(5, m) = B23
        end subroutine fill_shape

end program test_mpi_parcel_split
