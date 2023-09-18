! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose
    use constants, only : zero, two, one, f12, f13, f23
    use parcel_container, only : parcels, n_parcels, n_total_parcels, parcel_alloc
    use parcel_ellipsoid, only : get_abc, get_eigenvalues
    use parcel_split_mod, only : parcel_split
    use parcel_interpl, only : trilinear
    use parameters, only : dx, vcell, ncell,            &
                           extent, lower, nx, ny, nz,   &
                           max_num_parcels
    use timer, only : start_timer, stop_timer
    use field_mpi, only : field_halo_fill, field_halo_swap
    use dimensions, only : I_X, I_Y, I_Z
    use fields, only : vortg, tbuoyg
    use field_mpi, only : field_mpi_alloc                   &
                        , field_halo_fill                   &
                        , field_interior_to_buffer          &
                        , interior_to_halo_communication    &
                        , field_buffer_to_halo              &
                        , field_mpi_dealloc
#ifndef ENABLE_DRY_MODE
    use fields, only : humg
#endif
    use field_ops, only : get_mean, get_rms, get_abs_max
    use omp_lib
    use mpi_environment
    use mpi_layout, only : box
    use mpi_utils, only : mpi_print, mpi_exit_on_error
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    integer :: init_timer

    double precision, allocatable :: weights(:, :, :, :), apar(:)
    integer, allocatable :: is(:), js(:), ks(:)

    private :: weights, apar, is, js, ks

    double precision :: field_tol = 1.0e-12

    private :: init_refine,                 &
               alloc_and_precompute,        &
               dealloc

    contains

        ! This subroutine is only used in the unit test
        ! "test_parcel_init"
        subroutine unit_test_parcel_init_alloc
            call alloc_and_precompute
        end subroutine unit_test_parcel_init_alloc


        ! Allocate parcel container and sets values for parcel attributes
        ! to their default values.
        subroutine parcel_default
            double precision :: lam, l23
            integer          :: n

            call start_timer(init_timer)

            call parcel_alloc(max_num_parcels)

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel%n_per_cell * box%ncell

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                call mpi_exit_on_error
            endif

            n_total_parcels = n_parcels
            if (world%size > 1) then
                call mpi_blocking_reduce(n_total_parcels, MPI_SUM, world)
            endif

            call init_regular_positions

            ! initialize the volume of each parcel
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%volume(n) = vcell / dble(parcel%n_per_cell)
            enddo
            !$omp end do
            !$omp end parallel

            ! aspect ratio: lam = a / c
            lam = maxval((/dx(1) / dx(2), dx(2) / dx(1),   &
                           dx(1) / dx(3), dx(3) / dx(1),   &
                           dx(2) / dx(3), dx(3) / dx(2)/))

            !$omp parallel default(shared)
            !$omp do private(n, l23)
            do n = 1, n_parcels
                ! set all to zero
                parcels%B(:, n) = zero

                l23 = (lam * get_abc(parcels%volume(n))) ** f23

                ! B11
                parcels%B(1, n) = l23

                ! B22
                parcels%B(4, n) = l23
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(lam)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%vorticity(:, n) = zero
                parcels%buoyancy(n) = zero
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n) = zero
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(init_timer)

        end subroutine parcel_default


        ! Position parcels regularly in the domain.
        subroutine init_regular_positions
            integer          :: ix, i, iz, j, iy, k, l, n_per_dim
            double precision :: im, corner(3)

            ! number of parcels per dimension
            n_per_dim = int(dble(parcel%n_per_cell) ** f13)
            if (n_per_dim ** 3 .ne. parcel%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel%n_per_cell, ") not a cubic."
                stop
            endif

            im = one / dble(n_per_dim)

            l = 1
            do iz = 0, nz-1
                do iy = box%lo(2), box%hi(2)
                    do ix = box%lo(1), box%hi(1)
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do k = 1, n_per_dim
                            do j = 1, n_per_dim
                                do i = 1, n_per_dim
                                    parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                                    parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                                    parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                                    l = l + 1
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_parcels == l - 1) then
                print *, "Number of parcels disagree!"
                stop
            endif
        end subroutine init_regular_positions

        subroutine init_refine(lam)
            double precision, intent(inout) :: lam
            double precision                :: evals(3) ! = (a2, b2, c2)

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call parcel_split
                evals = get_eigenvalues(parcels%B(1, :), parcels%volume(1))
                lam = dsqrt(evals(1) / evals(3))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of trilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute
            double precision, allocatable :: resi(:, :, :)
            double precision              :: rsum
            integer                       :: n, i, j, k

            allocate(resi(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
            allocate(apar(n_parcels))
            allocate(weights(2, 2, 2, n_parcels))
            allocate(is(n_parcels))
            allocate(js(n_parcels))
            allocate(ks(n_parcels))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(l, n) reduction(+:resi)
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(:, n), is(n), js(n), ks(n), weights(:, :, :, n))

                i = is(n)
                j = js(n)
                k = ks(n)

                ! catch if in halo
                if ((k < 0) .or. (k >= nz)) then
                    print *, "Error: Tries to access undefined halo grid point."
                    stop
                endif

                resi(k:k+1, j:j+1, i:i+1) = resi(k:k+1, j:j+1, i:i+1) + weights(:, :, :, n)
            enddo
            !$omp end parallel do

            call field_halo_swap(resi)

            !Double edge values at iz = 0 and nz:
            resi(0,  :, :) = two * resi(0,  :, :)
            resi(nz, :, :) = two * resi(nz, :, :)

            ! Determine local inverse density of parcels (apar)
            !$omp parallel do default(shared) private(n, rsum)
            do n = 1, n_parcels
                rsum = sum(resi(ks(n):ks(n)+1, js(n):js(n)+1, is(n):is(n)+1) * weights(:, :, :, n))
                apar(n) = one / rsum
            enddo
            !$omp end parallel do

            deallocate(resi)

        end subroutine alloc_and_precompute

        subroutine dealloc
            deallocate(apar)
            deallocate(weights)
            deallocate(is)
            deallocate(js)
            deallocate(ks)
        end subroutine dealloc

        ! Initialise parcel attributes from gridded quantities.
        ! Attention: This subroutine currently only supports
        !            vorticity, buoyancy and humidity field.
        subroutine init_parcels_from_grids
            integer :: l

            call start_timer(init_timer)

            call alloc_and_precompute


            ! make sure halo grid points are filled
            call init_fill_halo

            do l = 1, 3
                call gen_parcel_scalar_attr(vortg(:, :, :, l), field_tol, parcels%vorticity(l, :))
            enddo

            call gen_parcel_scalar_attr(tbuoyg, field_tol, parcels%buoyancy)

#ifndef ENABLE_DRY_MODE
            call gen_parcel_scalar_attr(humg, field_tol, parcels%humidity)
#endif

            call dealloc

            call stop_timer(init_timer)

        end subroutine init_parcels_from_grids

        ! Generates the parcel attribute "par" from the field values provided
        ! in "field" (see Fontane & Dritschel, J. Comput. Phys. 2009, section 2.2)
        ! Precondition: The halo grid points of the field input must be filled.
        subroutine gen_parcel_scalar_attr(field, tol, par)
            double precision, intent(in)  :: field(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision, intent(in)  :: tol
            double precision, intent(out) :: par(:)
            double precision :: resi(-1:nz+1, box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1))
            double precision :: rms, rtol, rerr, rsum, fsum, avg_field
            integer          :: n, i, j, k

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    call mpi_print('Generate parcel attribute')
                endif
#endif

            ! Compute mean field value:
            avg_field = get_mean(field)

            resi(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) = &
                (field(0:nz, box%lo(2):box%hi(2), box%lo(1):box%hi(1)) - avg_field) ** 2

            ! "resi" is already squared, we only need to get the mean and apply the square root
            ! to get the rms
            rms = dsqrt(get_mean(resi))

            if (rms == zero) then
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_parcels
                    ! assign mean value
                    par(n) = avg_field
                enddo
                !$omp end do
                !$omp end parallel
                return
            endif

            ! Maximum error permitted below in gridded residue:
            rtol = rms * tol

            ! Initialise (volume-weighted) parcel attribute with a guess
            !$omp parallel do default(shared) private(n, fsum)
            do n = 1, n_parcels
                fsum = sum(field(ks(n):ks(n)+1, js(n):js(n)+1, is(n):is(n)+1) * weights(:, :, :, n))
                par(n) = apar(n) * fsum
            enddo
            !$omp end parallel do

            ! Iteratively compute a residual and update (volume-weighted) attribute:
            rerr = one

            do while (rerr .gt. rtol)
                !Compute residual:
                resi = zero
                do n = 1, n_parcels
                    i = is(n)
                    j = js(n)
                    k = ks(n)

                    resi(k:k+1, j:j+1, i:i+1) = resi(k:k+1, j:j+1, i:i+1) &
                                              + weights(:, :, :, n) * par(n)
                enddo

                call field_halo_swap(resi)

                resi(0, :, :)    = two * resi(0, :, :)
                resi(nz, :, :)   = two * resi(nz, :, :)
                resi(0:nz, :, :) = field(0:nz, :, :) - resi(0:nz, :, :)

                !Update (volume-weighted) attribute:
                !$omp parallel do default(shared) private(n, rsum, i, j, k)
                do n = 1, n_parcels
                    rsum = sum(resi(ks(n):ks(n)+1, js(n):js(n)+1, is(n):is(n)+1) * weights(:, :, :, n))
                    par(n) = par(n) + apar(n) * rsum
                enddo
                !$omp end parallel do

                !Compute maximum error:
                rerr = get_abs_max(resi)

#ifdef ENABLE_VERBOSE
                if (verbose .and. (world%rank == world%root)) then
                    print *, ' Max abs error = ', rerr
                endif
#endif
            enddo

            !Finally divide by parcel volume to define attribute:
            ! (multiply with vcell since algorithm is designed for volume fractions)
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                par(n) = vcell * par(n) / parcels%volume(n)
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine gen_parcel_scalar_attr

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_fill_halo
#ifndef ENABLE_DRY_MODE
            integer, parameter :: n_fields = 5
#else
            integer, parameter :: n_fields = 4
#endif
            call field_mpi_alloc(n_fields)

            call field_interior_to_buffer(vortg(:, :, :, I_X), 1)
            call field_interior_to_buffer(vortg(:, :, :, I_Y), 2)
            call field_interior_to_buffer(vortg(:, :, :, I_Z), 3)
            call field_interior_to_buffer(tbuoyg, 4)
#ifndef ENABLE_DRY_MODE
            call field_interior_to_buffer(humg, 5)
#endif

            call interior_to_halo_communication

            call field_buffer_to_halo(vortg(:, :, :, I_X), 1, .false.)
            call field_buffer_to_halo(vortg(:, :, :, I_Y), 2, .false.)
            call field_buffer_to_halo(vortg(:, :, :, I_Z), 3, .false.)
            call field_buffer_to_halo(tbuoyg, 4, .false.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_halo(humg, 5, .false.)
#endif
            call field_mpi_dealloc

        end subroutine init_fill_halo

end module parcel_init
