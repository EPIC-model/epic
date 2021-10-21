! =============================================================================
!                        Test parcel correction module
!
!         This unit test checks the parcel correction by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 500 relaxation steps.
! =============================================================================
program test_parcel_correction
    use unit_test
    use options, only : parcel
    use constants, only : pi, one, zero, f14, f32
    use parcel_container
    use parcel_correction
    use parcel_interpl, only : vol2grid
    use parcel_ellipse, only : get_ab
    use parcel_init, only : init_regular_positions
    use parameters, only : lower, extent, update_parameters, vcell, nx, ny, nz
    use fields, only : volg
    use timer

    implicit none

    double precision :: final_error, init_error, val, tmp
    integer :: i, n, sk, d
    integer, allocatable :: seed(:)
    double precision, parameter :: dev = 0.005d0

    call random_seed(size=sk)
    allocate(seed(1:sk))
    seed(:) = 42
    call random_seed(put=seed)

    call  parse_command_line

    call register_timer('laplace correction', lapl_corr_timer)
    call register_timer('gradient correction', grad_corr_timer)

    nx = 1
    ny = 32
    nz = 32
#ifdef ENABLE_3D
    nx = 32
    lower  = (/zero, zero, zero/)
    extent = (/one, one, one/)
#else
    lower  = (/zero, zero/)
    extent = (/one, one/)
#endif

    call update_parameters

    allocate(volg(-1:nz+1, 0:ny-1, 0:nx-1))

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)

    parcel%n_per_cell = 4
    call init_regular_positions

    ! add some deviation
    do n = 1, n_parcels
        do d = 1, ndim
            tmp = dev
            call random_number(val)
            if (val < 0.5) then
                tmp = -dev
            endif
            parcels%position(n, d) = parcels%position(n, d) + tmp
        enddo
    enddo

    volg = zero

    parcels%volume = f14 * vcell

#ifdef ENABLE_3D
    !FIXME
#else
    ! b11
    parcels%B(:, 1) = get_ab(parcels%volume(1:n_parcels))

    ! b12
    parcels%B(:, 2) = zero
#endif

    call vol2grid

    init_error = sum(abs(volg(0:nz, 0:ny-1, 0:nx-1) / vcell - one)) / (nx * ny * (nz+1))

    if (verbose) then
        write(*,*) 'Test laplace and gradient'
        write(*,*) 'iteration, average error, max absolute error'
        write(*,*) 0, init_error, maxval(abs(volg(0:nz, 0:ny-1, 0:nx-1) / vcell - one))
    endif

    call init_parcel_correction

    do i = 1, 20
        call apply_laplace
        call vol2grid
        call apply_gradient(1.80d0,0.5d0)
        if (verbose) then
            call vol2grid
            write(*,*) i, sum(abs(volg(0:nz, 0:ny-1, 0:nx-1) / vcell - one)) / (nx * ny * (nz+1)), &
                          maxval(abs(volg(0:nz, 0:ny-1, 0:nx-1) / vcell - one))
        endif
    enddo

    call vol2grid

    final_error = sum(abs(volg(0:nz, 0:ny-1, 0:nx-1) / vcell - one)) / (nx * ny * (nz+1))

    call print_result_dp('Test laplace and gradient', final_error, init_error)

    deallocate(volg)

end program test_parcel_correction
