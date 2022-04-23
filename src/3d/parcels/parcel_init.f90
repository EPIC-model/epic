! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12, f13, f23
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipsoid, only : get_abc, get_eigenvalues
    use parcel_split_mod, only : parcel_split
    use parcel_interpl, only : trilinear, ngp
    use parameters, only : dx, vcell, ncell,            &
                           extent, lower, nx, ny, nz,   &
                           max_num_parcels
    use timer, only : start_timer, stop_timer
    use field_mpi, only : field_halo_fill
    use omp_lib
    use mpi_layout, only : box
    implicit none

    integer :: init_timer

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :), ks(:, :)

    private :: weights, apar, is, js, ks


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
            double precision             :: lam, l23
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

            call init_from_grids(fname, tol)

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
                call parcel_split(parcels, parcel%lambda_max)
                evals = get_eigenvalues(parcels%B(1, :), parcels%volume(1))
                lam = dsqrt(evals(1) / evals(3))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of trilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute
            double precision, allocatable :: resi(:, :, :)
            double precision              :: rsum
            integer                       :: l, n

            allocate(resi(0:nz, 0:ny-1, 0:nx-1))
            allocate(apar(n_parcels))
            allocate(weights(ngp, n_parcels))
            allocate(is(ngp, n_parcels))
            allocate(js(ngp, n_parcels))
            allocate(ks(ngp, n_parcels))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(l, n) reduction(+:resi)
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(:, n), is(:, n), js(:, n), ks(:, n), weights(:, n))

                do l = 1, ngp
                    ! catch if in halo
                    if ((ks(l, n) < 0) .or. (ks(l, n) > nz)) then
                        print *, "Error: Tries to access undefined halo grid point."
                        stop
                    endif
                    resi(ks(l, n), js(l, n), is(l, n)) = resi(ks(l, n), js(l, n), is(l, n)) + weights(l, n)
                enddo
            enddo
            !$omp end parallel do

            !Double edge values at iz = 0 and nz:
            resi(0,  :, :) = two * resi(0,  :, :)
            resi(nz, :, :) = two * resi(nz, :, :)

            ! Determine local inverse density of parcels (apar)
            !$omp parallel do default(shared) private(l, n, rsum)
            do n = 1, n_parcels
                rsum = zero
                do l = 1, ngp
                    rsum = rsum + resi(ks(l, n), js(l, n), is(l, n)) * weights(l, n)
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
            deallocate(ks)
        end subroutine dealloc

        ! Initialise parcel attributes from gridded quantities.
        ! Attention: This subroutine currently only supports
        !            vorticity and buoyancy fields.
        subroutine init_from_grids(ncfname, tol)
            use netcdf_reader
            character(*),     intent(in)  :: ncfname
            double precision, intent(in)  :: tol
            double precision, allocatable :: buffer(:, :, :)
            integer                       :: ncid
            integer                       :: n_steps, start(4), cnt(4)
            integer                       :: lo(3), hi(3)

            call alloc_and_precompute

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)


            ! allocate with halo grid points
            allocate(buffer(box%hlo(3):box%hhi(3), &
                            box%hlo(2):box%hhi(2), &
                            box%hlo(1):box%hhi(1)))

            ! read without halo grid points
            ! we must add +1 since index starts at 1
            lo = box%lo
            hi = box%hi
            start(1:3) = lo + 1
            start(4)   = n_steps

            cnt(1:3) = hi - lo + 1
            cnt(4)   = 1

            if (has_dataset(ncid, 'x_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'x_vorticity', buffer(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                         start=start, cnt=cnt)

                call field_halo_fill(buffer)

                call gen_parcel_scalar_attr(buffer, tol, parcels%vorticity(1, :))
            endif

            if (has_dataset(ncid, 'y_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'y_vorticity', buffer(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                         start=start, cnt=cnt)

                call field_halo_fill(buffer)

                call gen_parcel_scalar_attr(buffer, tol, parcels%vorticity(2, :))
            endif

            if (has_dataset(ncid, 'z_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'z_vorticity', buffer(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                         start=start, cnt=cnt)

                call field_halo_fill(buffer)

                call gen_parcel_scalar_attr(buffer, tol, parcels%vorticity(3, :))
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'buoyancy', buffer(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                         start=start, cnt=cnt)

                call field_halo_fill(buffer)

                call gen_parcel_scalar_attr(buffer, tol, parcels%buoyancy)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'humidity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'humidity', buffer(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)), &
                                         start=start, cnt=cnt)

                call field_halo_fill(buffer)

                call gen_parcel_scalar_attr(buffer, tol, parcels%humidity)
            endif
#endif
            call close_netcdf_file(ncid)

            deallocate(buffer)

            call dealloc

        end subroutine init_from_grids

        ! Generates the parcel attribute "par" from the field values provided
        ! in "field" (see Fontane & Dritschel, J. Comput. Phys. 2009, section 2.2)
        subroutine gen_parcel_scalar_attr(field, tol, par)
            double precision, intent(in)  :: field(-1:nz+1, 0:ny-1, 0:nx-1)
            double precision, intent(in)  :: tol
            double precision, intent(out) :: par(:)
            double precision :: resi(0:nz, 0:ny-1, 0:nx-1)
            double precision :: rms, rtol, rerr, rsum, fsum, avg_field
            integer          :: l, n

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
            !$omp parallel do default(shared) private(l, n, fsum)
            do n = 1, n_parcels
                fsum = zero
                do l = 1, ngp
                    fsum = fsum + field(ks(l, n), js(l, n), is(l, n)) * weights(l, n)
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
                        resi(ks(l, n), js(l, n), is(l, n)) = resi(ks(l, n), js(l, n), is(l, n)) &
                                                           + weights(l, n) * par(n)
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
                        rsum = rsum + resi(ks(l, n), js(l, n), is(l, n)) * weights(l, n)
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
