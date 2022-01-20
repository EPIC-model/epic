! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12, max_num_parcels
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parcel_interpl, only : bilinear, ngp
    use parameters, only : update_parameters,   &
                           dx, vcell, ncell,    &
                           extent, lower, nx, nz
    use h5_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: init_timer

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :)

    private :: weights, apar, is, js


    private :: init_refine,                 &
               init_from_grids,             &
               fill_field_from_buffer_2d,   &
               alloc_and_precompute,        &
               dealloc

    contains

        ! This subroutine is only used in the unit test
        ! "test_parcel_init"
        subroutine unit_test_parcel_init_alloc
            call alloc_and_precompute
        end subroutine unit_test_parcel_init_alloc

        ! This subroutine reads parcels from an EPIC output file
        subroutine read_parcels(h5fname, step)
            character(*),     intent(in)  :: h5fname
            integer,          intent(in)  :: step
            integer(hid_t)                :: h5handle, group
            integer                       :: ncells(2)
            character(:), allocatable     :: grn
            double precision, allocatable :: buffer_1d(:),   &
                                             buffer_2d(:, :)
            logical                       :: l_valid = .false.


            call start_timer(init_timer)

            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)

            ! read domain dimensions
            call read_h5_box(h5handle, ncells, extent, lower)
            nx = ncells(1)
            nz = ncells(2)

            ! update global parameters
            call update_parameters

            grn = trim(get_step_group_name(step))
            call open_h5_group(h5handle, grn, group)
            call get_num_parcels(group, n_parcels)

            if (n_parcels > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_parcels, ". Exiting."
                stop
            endif

            ! Be aware that the starting index of buffer_1d and buffer_2d
            ! is 0; hence, the range is 0:n_parcels-1 in contrast to the
            ! parcel container where it is 1:n_parcels.

            if (has_dataset(group, 'B')) then
                call read_h5_dataset(group, 'B', buffer_2d)
                parcels%B(:, 1:n_parcels) = buffer_2d
                deallocate(buffer_2d)
            else
                print *, "The parcel shape must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'position')) then
                call read_h5_dataset(group, 'position', buffer_2d)
                parcels%position(:, 1:n_parcels) = buffer_2d
                deallocate(buffer_2d)
            else
                print *, "The parcel position must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'volume')) then
                call read_h5_dataset(group, 'volume', buffer_1d)
                parcels%volume(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            else
                print *, "The parcel volume must be present! Exiting."
                stop
            endif

            if (has_dataset(group, 'vorticity')) then
                l_valid = .true.
                call read_h5_dataset(group, 'vorticity', buffer_1d)
                parcels%vorticity(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            endif

            if (has_dataset(group, 'buoyancy')) then
                l_valid = .true.
                call read_h5_dataset(group, 'buoyancy', buffer_1d)
                parcels%buoyancy(1:n_parcels) = buffer_1d
                deallocate(buffer_1d)
            endif

            if (.not. l_valid) then
                print *, "Either the parcel buoyancy or vorticity must be present! Exiting."
                stop
            endif

            call close_h5_group(group)
            call close_h5_file(h5handle)

            call stop_timer(init_timer)

        end subroutine read_parcels

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(h5fname, tol)
            character(*),     intent(in) :: h5fname
            double precision, intent(in) :: tol
            double precision             :: lam, ratio
            integer(hid_t)               :: h5handle
            integer                      :: n, ncells(2)

            call start_timer(init_timer)

            ! read domain dimensions
            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)
            call read_h5_box(h5handle, ncells, extent, lower)
            nx = ncells(1)
            nz = ncells(2)
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

            ratio = dx(1) / dx(2)

            ! aspect ratio: lam = a / b
            lam = max(dx(2) / dx(1), ratio)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                ! B11
                parcels%B(n, 1) = ratio * get_ab(parcels%volume(n))

                ! B12
                parcels%B(n, 2) = zero
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(lam)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                parcels%vorticity(n) = zero
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
            integer          :: ix, i, iz, j, k, n_per_dim
            double precision :: im, corner(2)

            ! number of parcels per dimension
            n_per_dim = int(dsqrt(dble(parcel%n_per_cell)))
            if (n_per_dim ** 2 .ne. parcel%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel%n_per_cell, ") not a square."
                stop
            endif

            im = one / dble(n_per_dim)

            k = 1
            do iz = 0, nz-1
                do ix = 0, nx-1
                    corner = lower + dble((/ix, iz/)) * dx
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            parcels%position(k, 1) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(k, 2) = corner(2) + dx(2) * (dble(j) - f12) * im
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_parcels == k - 1) then
                print *, "Number of parcels disagree!"
                stop
            endif
        end subroutine init_regular_positions

        subroutine init_refine(lam)
            double precision, intent(inout) :: lam
            double precision                :: B22, a2

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call split_ellipses(parcels, parcel%lambda_max)
                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of bilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute
            double precision :: resi(0:nz, 0:nx-1), rsum
            integer          :: n, l

            allocate(apar(n_parcels))
            allocate(weights(n_parcels, ngp))
            allocate(is(n_parcels, ngp))
            allocate(js(n_parcels, ngp))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(n, l) reduction(+:resi)
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call bilinear(parcels%position(n, :), is(n, :), js(n, :), weights(n, :))

                do l = 1, ngp
                    resi(js(n, l), is(n, l)) = resi(js(n, l), is(n, l)) + weights(n, l)
                enddo
            enddo
            !$omp end parallel do

            !Double edge values at iz = 0 and nz:
            resi(0, :) = two * resi(0, :)
            resi(nz,:) = two * resi(nz,:)

            ! Determine local inverse density of parcels (apar)
            !$omp parallel do default(shared) private(n, l, rsum)
            do n = 1, n_parcels
                rsum = zero
                do l = 1, ngp
                    rsum = rsum + resi(js(n, l), is(n, l)) * weights(n, l)
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
        end subroutine dealloc


        ! Initialise parcel attributes from gridded quantities.
        ! Attention: This subroutine currently only supports
        !            vorticity and buoyancy fields.
        subroutine init_from_grids(h5fname, tol)
            character(*),     intent(in)  :: h5fname
            double precision, intent(in)  :: tol
            double precision, allocatable :: buffer_2d(:, :)
            double precision              :: field_2d(-1:nz+1, 0:nx-1)
            integer(hid_t)                :: h5handle, group
            integer                       :: n_steps
            character(:), allocatable     :: grn

            call alloc_and_precompute

            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)

            call get_num_steps(h5handle, n_steps)

            group = h5handle
            if (n_steps > 0) then
                ! we need to subtract 1 since it starts with 0
                grn = trim(get_step_group_name(n_steps - 1))
                call open_h5_group(h5handle, grn, group)
            endif

            if (has_dataset(group, 'vorticity')) then
                call read_h5_dataset(group, 'vorticity', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, tol, parcels%vorticity)
            endif


            if (has_dataset(group, 'buoyancy')) then
                call read_h5_dataset(group, 'buoyancy', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, tol, parcels%buoyancy)
            endif

            if (n_steps > 0) then
                call close_h5_group(group)
            endif
            call close_h5_file(h5handle)

            call dealloc

        end subroutine init_from_grids

        ! After reading the H5 dataset into the buffer, copy
        ! the data to a field container
        ! @pre field and buffer must be of rank 2
        subroutine fill_field_from_buffer_2d(buffer, field)
            double precision, allocatable :: buffer(:, :)
            double precision              :: field(-1:nz+1, 0:nx-1)
            integer                       :: dims(2), bdims(2), i, j

            dims = (/nz+1, nx/)

            bdims = shape(buffer)
            if (.not. sum(dims - bdims) == 0) then
                print "(a32, i4, a1, i4, a6, i4, a1, i4, a1)", &
                      "Field dimensions do not agree: (", dims(1), ",", &
                      dims(2), ") != (", bdims(1), ",", bdims(2), ")"
                stop
            endif

            do j = 0, nz
                do i = 0, nx-1
                    field(j, i) = buffer(j, i)
                enddo
            enddo
        end subroutine fill_field_from_buffer_2d


        ! Generates the parcel attribute "par" from the field values provided
        ! in "field" (see Fontane & Dritschel, J. Comput. Phys. 2009, section 2.2)
        subroutine gen_parcel_scalar_attr(field, tol, par)
            double precision, intent(in)  :: field(-1:nz+1, 0:nx-1)
            double precision, intent(in)  :: tol
            double precision, intent(out) :: par(:)
            double precision :: resi(0:nz, 0:nx-1)
            double precision :: rms, rtol, rerr, rsum, fsum, avg_field
            integer          :: n, l

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, 'Generate parcel attribute'
                endif
#endif

            ! Compute mean field value:
            ! (divide by ncell since lower and upper edge weights are halved)
            avg_field = (f12 * sum(field(0, :) + field(nz, :)) &
                             + sum(field(1:nz-1,:))) / dble(ncell)

            resi(0:nz,:) = (field(0:nz,:) - avg_field) ** 2

            rms = dsqrt((f12 * sum(resi(0, :) + resi(nz, :)) &
                             + sum(resi(1:nz-1,:))) / dble(ncell))


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
                    fsum = fsum + field(js(n, l), is(n, l)) * weights(n, l)
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
                        resi(js(n, l), is(n, l)) = resi(js(n, l), is(n, l)) &
                                                 + weights(n, l) * par(n)
                    enddo
                enddo

                resi(0, :)    = two * resi(0, :)
                resi(nz, :)   = two * resi(nz, :)
                resi(0:nz, :) = field(0:nz, :) - resi(0:nz, :)

                !Update (volume-weighted) attribute:
                !$omp parallel do default(shared) private(n, rsum, l)
                do n = 1, n_parcels
                    rsum = zero
                    do l = 1, ngp
                        rsum = rsum + resi(js(n, l), is(n, l)) * weights(n, l)
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
