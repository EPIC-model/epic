! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_eigenvalue, get_B22
    use parcel_split, only : split_ellipses
    use parcel_interpl, only : bilinear, ngp
    use parameters, only : dx, vcell, ncell,        &
                           extent, lower, nx, nz,   &
                           max_num_parcels
    use netcdf_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: init_timer

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :)

    private :: weights, apar, is, js


    private :: init_refine,                 &
               init_from_grids,             &
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
        subroutine init_parcels(fname, tol)
            character(*),     intent(in) :: fname
            double precision, intent(in) :: tol
            double precision             :: lam, ratio
            integer                      :: n

            call start_timer(init_timer)

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
                parcels%B(1, n) = ratio * get_ab(parcels%volume(n))

                ! B12
                parcels%B(2, n) = zero

                ! B22
                parcels%B(3, n) = get_B22(parcels%B(1, n), zero, parcels%volume(1))
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

            call init_from_grids(fname, tol)

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
                            parcels%position(1, k) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(2, k) = corner(2) + dx(2) * (dble(j) - f12) * im
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
            double precision                :: a2

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call split_ellipses(parcels, parcel%lambda_max)
                a2 = get_eigenvalue(parcels%B(1, 1), parcels%B(2, 1), parcels%B(3, 1))
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of bilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute
            double precision, allocatable :: resi(:, :)
            double precision :: rsum
            integer          :: n, l

            allocate(apar(n_parcels))
            allocate(weights(ngp, n_parcels))
            allocate(is(ngp, n_parcels))
            allocate(js(ngp, n_parcels))
            allocate(resi(0:nz, 0:nx-1))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(n, l) reduction(+:resi)
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call bilinear(parcels%position(:, n), is(:, n), js(:, n), weights(:, n))

                do l = 1, ngp
                    ! catch if in halo
                    if ((js(l, n) < 0) .or. (js(l, n) > nz)) then
                        print *, "Error: Tries to access undefined halo grid point."
                        stop
                    endif
                    resi(js(l, n), is(l, n)) = resi(js(l, n), is(l, n)) + weights(l, n)
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
                    rsum = rsum + resi(js(l, n), is(l, n)) * weights(l, n)
                enddo
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
        end subroutine dealloc


        ! Initialise parcel attributes from gridded quantities.
        ! Attention: This subroutine currently only supports
        !            vorticity and buoyancy fields.
        subroutine init_from_grids(ncfname, tol)
            character(*),     intent(in)  :: ncfname
            double precision, intent(in)  :: tol
            double precision              :: buffer(-1:nz+1, 0:nx-1)
            integer                       :: ncid
            integer                       :: n_steps, start(3), cnt(3)

            call alloc_and_precompute

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            cnt  =  (/ nx, nz+1, 1       /)
            start = (/ 1,  1,    n_steps /)


            if (has_dataset(ncid, 'vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'vorticity', buffer(0:nz, :), start=start, cnt=cnt)
                call gen_parcel_scalar_attr(buffer, tol, parcels%vorticity)
            endif


            if (has_dataset(ncid, 'buoyancy')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'buoyancy', buffer(0:nz, :), start=start, cnt=cnt)
                call gen_parcel_scalar_attr(buffer, tol, parcels%buoyancy)
            endif

            call close_netcdf_file(ncid)

            call dealloc

        end subroutine init_from_grids

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
                    fsum = fsum + field(js(l, n), is(l, n)) * weights(l, n)
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
                        resi(js(l, n), is(l, n)) = resi(js(l, n), is(l, n)) &
                                                 + weights(l, n) * par(n)
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
                        rsum = rsum + resi(js(l, n), is(l, n)) * weights(l, n)
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
