! =============================================================================
!                       Test diverge module
!
!         This unit test checks the diverge module by initializing
!         the parcels with a small deviation from the optimal position.
!         It then performs 500 relaxation steps.
! =============================================================================
program test_diverge_gradient
    use unit_test
    use constants, only : pi
    use parcel_container
    use parcel_diverge
    use parcel_interpl, only : par2grid, grid2par
    use options, only : parcel_info, grid, interpl
    use ellipse, only : get_ab
    use parameters, only : extent, lower, update_parameters, vcell, dx, nx, nz, ngrid
    implicit none

    double precision :: volg(-1:33, 0:31, 1), final_error, init_error
    integer :: i, j, k, jj, ii

    call parse_command_line

    grid = (/33, 33/)

    extent =  (/0.4, 0.4/)

    call update_parameters()

    n_parcels = 4*nx*nz
    call parcel_alloc(n_parcels)


    k = 1
    do j = 0, nz-1
        do i = 0, nx-1
            do jj = 1, 4, 2
                do ii = 1, 4, 2
                    parcels%position(k, 1) = lower(1) + i * dx(1) + 0.25 * dx(1) * ii
                    parcels%position(k, 2) = lower(2) + j * dx(2) + 0.25 * dx(2) * jj

                    ! add some deviation
                    parcels%position(k, 1) = (1.0 + 1.0e-7) * parcels%position(k, 1)
                    parcels%position(k, 2) = (1.0 + 1.0e-7) * parcels%position(k, 2)
                    k = k + 1
                enddo
            enddo
        enddo
    enddo

    volg = 0.0

    parcel_info%is_elliptic = .true.

    parcels%volume = 0.25d0 * vcell

    ! b11
    parcels%B(:, 1) = get_ab(parcels%volume(1:n_parcels, 1))

    ! b12
    parcels%B(:, 2) = 0.0d0

    call par2grid(parcels, parcels%volume, volg)

    init_error = sum(abs(volg(0:nz, 0:nx-1, :) - vcell))

    if (verbose) then
        write(*,*) 'Test diverge and gradient, initial error'
        write(*,*) init_error
    endif

    call init_diverge

    do i = 1, 500
        call par2grid(parcels, parcels%volume, volg)
        call apply_diverge(volg)
        call par2grid(parcels, parcels%volume, volg)
        call apply_gradient(volg,0.30d0)
        call par2grid(parcels, parcels%volume, volg)
        if (verbose) then
            write(*,*) 'Test diverge and gradient, error after iteration ', i
            write(*,*) sum(abs(volg(0:nz, 0:nx-1, :) - vcell))
        endif
    enddo

    if (verbose) then
        write(*,*) 'Test diverge and gradient, final error'
        write(*,*) final_error
    endif

    call print_result_dp('Test diverge and gradient', final_error, init_error)

end program test_diverge_gradient
