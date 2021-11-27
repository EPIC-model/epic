! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12, f13, f23, max_num_parcels
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid, only : get_abc, get_eigenvalues
    use parcel_split_mod, only : parcel_split
    use parcel_interpl, only : trilinear, ngp
    use parameters, only : update_parameters,   &
                           dx, vcell, ncell,    &
                           extent, lower, nx, ny, nz
    use h5_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: init_timer

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :), ks(:, :)

    private :: weights, apar, is, js, ks


    private :: init_refine,                 &
               init_from_grids,             &
               fill_field_from_buffer_3d,   &
               fill_field_from_buffer_4d,   &
               alloc_and_precompute,        &
               dealloc

    contains

        ! This subroutine is only used in the unit test
        ! "test_parcel_init"
        subroutine unit_test_parcel_init_alloc
            call alloc_and_precompute
        end subroutine unit_test_parcel_init_alloc


        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(h5fname, tol)
            character(*),     intent(in) :: h5fname
            double precision, intent(in) :: tol
            double precision             :: lam, l23
            integer(hid_t)               :: h5handle
            integer                      :: n, ncells(3)

            call start_timer(init_timer)

            ! read domain dimensions
            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)
            call read_h5_box(h5handle, ncells, extent, lower)
            nx = ncells(1)
            ny = ncells(2)
            nz = ncells(3)
            call close_h5_file(h5handle)

            ! update global parameters
            call update_parameters

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel%n_per_cell * ncell

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                stop
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
                parcels%B(n, :) = zero

                l23 = (lam * get_abc(parcels%volume(n))) ** f23

                ! B11
                parcels%B(n, 1) = l23

                ! B22
                parcels%B(n, 4) = l23
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(lam)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%vorticity(n, :) = zero
                parcels%buoyancy(n) = zero
#ifndef ENABLE_DRY_MODE
                parcels%humidity(n) = zero
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call init_from_grids(h5fname, tol)

            call stop_timer(init_timer)

        end subroutine init_parcels


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
                do iy = 0, ny-1
                    do ix = 0, nx-1
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do k = 1, n_per_dim
                            do j = 1, n_per_dim
                                do i = 1, n_per_dim
                                    parcels%position(l, 1) = corner(1) + dx(1) * (dble(i) - f12) * im
                                    parcels%position(l, 2) = corner(2) + dx(2) * (dble(j) - f12) * im
                                    parcels%position(l, 3) = corner(3) + dx(3) * (dble(k) - f12) * im
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
                call parcel_split(parcels, parcel%lambda_max)
                evals = get_eigenvalues(parcels%B(1, :), parcels%volume(1))
                lam = dsqrt(evals(1) / evals(3))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of trilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute
            double precision :: resi(0:nz, 0:ny-1, 0:nx-1), rsum
            integer          :: n, l

            allocate(apar(n_parcels))
            allocate(weights(n_parcels, ngp))
            allocate(is(n_parcels, ngp))
            allocate(js(n_parcels, ngp))
            allocate(ks(n_parcels, ngp))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(n, l) reduction(+:resi)
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(n, :), is(n, :), js(n, :), ks(n, :), weights(n, :))

                do l = 1, ngp
                    resi(ks(n, l), js(n, l), is(n, l)) = resi(ks(n, l), js(n, l), is(n, l)) + weights(n, l)
                enddo
            enddo
            !$omp end parallel do

            !Double edge values at iz = 0 and nz:
            resi(0,  :, :) = two * resi(0,  :, :)
            resi(nz, :, :) = two * resi(nz, :, :)

            ! Determine local inverse density of parcels (apar)
            !$omp parallel do default(shared) private(n, l, rsum)
            do n = 1, n_parcels
                rsum = zero
                do l = 1, ngp
                    rsum = rsum + resi(ks(n, l), js(n, l), is(n, l)) * weights(n, l)
                enddo
                apar(n) = one / rsum
            enddo
            !$omp end parallel do

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
        !            vorticity and buoyancy fields.
        subroutine init_from_grids(h5fname, tol)
            character(*),     intent(in)  :: h5fname
            double precision, intent(in)  :: tol
            double precision, allocatable :: buffer_3d(:, :, :), buffer_4d(:, :, :, :)
            double precision              :: field_3d(-1:nz+1, 0:ny-1, 0:nx-1)
            integer(hid_t)                :: h5handle
            integer                       :: l

            call alloc_and_precompute

            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)

            if (has_dataset(h5handle, 'vorticity')) then
                call read_h5_dataset(h5handle, 'vorticity', buffer_4d)
                call fill_field_from_buffer_4d(buffer_4d, buffer_4d)
                deallocate(buffer_4d)
                do l = 1, 3
                    call gen_parcel_scalar_attr(buffer_4d(:, :, :, l), tol, parcels%vorticity(:, l))
                enddo
            endif


            if (has_dataset(h5handle, 'buoyancy')) then
                call read_h5_dataset(h5handle, 'buoyancy', buffer_3d)
                call fill_field_from_buffer_3d(buffer_3d, field_3d)
                deallocate(buffer_3d)
                call gen_parcel_scalar_attr(field_3d, tol, parcels%buoyancy)
            endif

            call close_h5_file(h5handle)

            call dealloc

        end subroutine init_from_grids

        ! After reading the H5 dataset scalar field into the buffer, copy
        ! the data to a field container
        ! @pre field and buffer must be of rank 3
        subroutine fill_field_from_buffer_3d(buffer, field)
            double precision, allocatable :: buffer(:, :, :)
            double precision              :: field(-1:nz+1, 0:ny-1, 0:nx-1)
            integer                       :: dims(3), bdims(3), i, j, k

            dims = (/nz+1, ny, nx/)

            bdims = shape(buffer)
            if (.not. sum(dims - bdims) == 0) then
                print "(a32, i4, a1, i4, a1, i4, a6, i4, a1, i4, a1, i4, a1)", &
                      "Field dimensions do not agree: (", dims(1), ",", &
                      dims(2), ",", dims(3), ") != (", bdims(1), ",", bdims(2), ",", bdims(3), ")"
                stop
            endif

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz
                        field(k, j, i) = buffer(k, j, i)
                    enddo
                enddo
            enddo
        end subroutine fill_field_from_buffer_3d

        ! After reading the H5 dataset vector field into the buffer, copy
        ! the data to a field container
        ! @pre field and buffer must be of rank 4
        subroutine fill_field_from_buffer_4d(buffer, field)
            double precision, allocatable :: buffer(:, :, :, :)
            double precision              :: field(-1:nz+1, 0:ny-1, 0:nx-1, 3)
            integer                       :: dims(4), bdims(4), i, j, k

            dims = (/nz+1, ny, nx, 3/)

            bdims = shape(buffer)
            if (.not. sum(dims - bdims) == 0) then
                print "(a32, i4, a1, i4, a1, i4, a1, i4, a6, i4, a1, i4, a1, i4, a1, i4, a1)", &
                      "Field dimensions do not agree: (", dims(1), ",", &
                      dims(2), ",", dims(3), ",", dims(4), ") != (",    &
                      bdims(1), ",", bdims(2), ",", bdims(3), ",", bdims(4), ")"
                stop
            endif

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz
                        field(k, j, i, :) = buffer(k, j, i, :)
                    enddo
                enddo
            enddo
        end subroutine fill_field_from_buffer_4d

        ! Generates the parcel attribute "par" from the field values provided
        ! in "field" (see Fontane & Dritschel, J. Comput. Phys. 2009, section 2.2)
        subroutine gen_parcel_scalar_attr(field, tol, par)
            double precision, intent(in)  :: field(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision, intent(in)  :: tol
            double precision, intent(out) :: par(:)
            double precision :: resi(0:nz, 0:ny-1, 0:nx-1)
            double precision :: rms, rtol, rerr, rsum, fsum, avg_field
            integer          :: n, l

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, 'Generate parcel attribute'
                endif
#endif

            ! Compute mean field value:
            ! (divide by ncell since lower and upper edge weights are halved)
            avg_field = (f12 * sum(field(0, :, :) + field(nz, :, :)) &
                             + sum(field(1:nz-1, :, :))) / dble(ncell)

            resi(0:nz, :, :) = (field(0:nz, :, :) - avg_field) ** 2

            rms = dsqrt((f12 * sum(resi(0, :, :) + resi(nz, :, :)) &
                             + sum(resi(1:nz-1, :, :))) / dble(ncell))


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
            !$omp parallel do default(shared) private(n, l, fsum)
            do n = 1, n_parcels
                fsum = zero
                do l = 1, ngp
                    fsum = fsum + field(ks(n, l), js(n, l), is(n, l)) * weights(n, l)
                enddo
                par(n) = apar(n) * fsum
            enddo
            !$omp end parallel do

            ! Iteratively compute a residual and update (volume-weighted) attribute:
            rerr = one

            do while (rerr .gt. rtol)
                !Compute residual:
                resi = zero
                do n = 1, n_parcels
                    do l = 1, ngp
                        resi(ks(n, l), js(n, l), is(n, l)) = resi(ks(n, l), js(n, l), is(n, l)) &
                                                           + weights(n, l) * par(n)
                    enddo
                enddo

                resi(0, :, :)    = two * resi(0, :, :)
                resi(nz, :, :)   = two * resi(nz, :, :)
                resi(0:nz, :, :) = field(0:nz, :, :) - resi(0:nz, :, :)

                !Update (volume-weighted) attribute:
                !$omp parallel do default(shared) private(n, rsum, l)
                do n = 1, n_parcels
                    rsum = zero
                    do l = 1, ngp
                        rsum = rsum + resi(ks(n, l), js(n, l), is(n, l)) * weights(n, l)
                    enddo
                    par(n) = par(n) + apar(n) * rsum
                enddo
                !$omp end parallel do

                !Compute maximum error:
                rerr = maxval(dabs(resi))

#ifdef ENABLE_VERBOSE
                if (verbose) then
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

end module parcel_init
