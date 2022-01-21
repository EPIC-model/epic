! =============================================================================
!                       Test parcel laplace correction module
!
!         This unit test checks the correction module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 20 relaxation steps.
! =============================================================================
program test_laplace_correction
    use unit_test
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f32
    use parcel_container
    use parcel_correction
    use parcel_interpl, only : vol2grid
    use parcel_ellipse, only : get_ab
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, nz
    use fields, only : volg
    use timer
    implicit none

    double precision :: final_error, init_error, val, tmp
    integer :: i, n, sk
    integer, allocatable :: seed(:)
    double precision, parameter :: dev = 0.005d0

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)

    call  parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)


    nx = 32
    nz = 32
    lower  = (/zero, zero/)
    extent = (/one, one/)

    call update_parameters

    allocate(volg(-1:nz+1, 0:nx-1))

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)

    parcel%n_per_cell = 4
    call init_regular_positions

    ! add some deviation
    do n = 1, n_parcels
        tmp = dev
        call random_number(val)
        if (val < 0.5) then
            tmp = -dev
        endif
        parcels%position(1, n) = parcels%position(1, n) + tmp

        tmp = dev
        call random_number(val)
        if (val < 0.5) then
            tmp = -dev
        endif
        parcels%position(2, n) = parcels%position(2, n) + tmp
    enddo

    volg = zero

    parcels%volume = f14 * vcell

    ! b11
    parcels%B(1, :) = get_ab(parcels%volume(1:n_parcels))

    ! b12
    parcels%B(2, :) = zero

    call vol2grid

    init_error = sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1))

    if (verbose) then
        write(*,*) 'test laplace correction'
        write(*,*) 'iteration, average error, max absolute error'
        write(*,*) 0, init_error, maxval(abs(volg(0:nz, 0:nx-1) / vcell - one))
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        if (verbose) then
            call vol2grid
            write(*,*) i, sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1)), &
                          maxval(abs(volg(0:nz, 0:nx-1) / vcell - one))
        endif
    enddo

    call vol2grid

    final_error = sum(abs(volg(0:nz, 0:nx-1) / vcell - one)) / (nx * (nz+1))

    call print_result_dp('Test laplace correction', final_error, init_error)

    deallocate(volg)

end program test_laplace_correction

