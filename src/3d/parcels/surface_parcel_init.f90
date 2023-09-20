! =============================================================================
!             This module initializes surface parcel default values.
! =============================================================================
module surface_parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, one, f12
    use surface_parcel_container, only : surface_parcel_container_type      &
                                       , lo_surf_parcels, n_lo_surf_parcels &
                                       , up_surf_parcels, n_up_surf_parcels
    use parcel_ellipse, only : get_ab, get_eigenvalue
    use surface_parcel_split, only : do_ellipse_split
    use surface_parcel_interpl, only : bilinear, ngp
    use parameters, only : dx, acell,                   &
                           extent, lower, nx, ny, nz,   &
                           max_num_surf_parcels
    use netcdf_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: surf_init_timer

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :)

    private :: weights, apar, is, js


    private :: init_refine,                 &
               init_parcels,                &
               init_from_grids,             &
               alloc_and_precompute,        &
               gen_parcel_scalar_attr,      &
               dealloc

    contains

!         ! This subroutine is only used in the unit test
!         ! "test_parcel_init"
!         subroutine unit_test_parcel_init_alloc
!             call alloc_and_precompute
!         end subroutine unit_test_parcel_init_alloc

        subroutine init_surface_parcels(fname, tol)
            character(*),     intent(in) :: fname
            double precision, intent(in) :: tol
            call init_parcels(lo_surf_parcels, n_lo_surf_parcels, 'lo', fname, tol)
            call init_parcels(up_surf_parcels, n_up_surf_parcels, 'up', fname, tol)
        end subroutine init_surface_parcels

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(s_parcels, n_par, which, fname, tol)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            character(2),                        intent(in)    :: which
            character(*),     intent(in)                       :: fname
            double precision, intent(in)                       :: tol
            double precision                                   :: lam, ratio
            integer                                            :: n

            call start_timer(surf_init_timer)

            ! set the number of parcels (see parcels.f90)
            ! we use "n_surf_per_cell" parcels per grid cell
            n_par = parcel%n_surf_per_cell * nx * ny

            if (n_par > max_num_surf_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_surf_parcels, ". Exiting."
                stop
            endif


            call init_regular_positions(s_parcels, n_par)

            ! initialize the area of each parcel
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                s_parcels%area(n) = acell / dble(parcel%n_surf_per_cell)
            enddo
            !$omp end do
            !$omp end parallel

            ratio = dx(1) / dx(2)

            ! aspect ratio: lam = a / b
            lam = max(dx(2) / dx(1), ratio)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                ! B11
                s_parcels%B(1, n) = ratio * get_ab(s_parcels%area(n))

                ! B12
                s_parcels%B(2, n) = zero

                ! B22
                s_parcels%B(3, n) = get_ab(s_parcels%area(n)) / ratio
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(s_parcels, n_par, lam)

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                s_parcels%buoyancy(n) = zero
                s_parcels%vorticity(:, n) = zero
#ifndef ENABLE_DRY_MODE
                s_parcels%humidity(n) = zero
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call init_from_grids(s_parcels, n_par, which, fname, tol)

            call stop_timer(surf_init_timer)

        end subroutine init_parcels


        ! Position parcels regularly in the domain.
        subroutine init_regular_positions(s_parcels, n_par)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            integer                                            :: ix, i, iy, j, k, n_per_dim
            double precision                                   :: im, corner(2)

            ! number of parcels per dimension
            n_per_dim = int(dsqrt(dble(parcel%n_surf_per_cell)))
            if (n_per_dim ** 2 .ne. parcel%n_surf_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel%n_surf_per_cell, ") not a square."
                stop
            endif

            im = one / dble(n_per_dim)

            k = 1
            do iy = 0, ny-1
                do ix = 0, nx-1
                    corner = lower(1:2) + dble((/ix, iy/)) * dx(1:2)
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            s_parcels%position(1, k) = corner(1) + dx(1) * (dble(i) - f12) * im
                            s_parcels%position(2, k) = corner(2) + dx(2) * (dble(j) - f12) * im
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_par == k - 1) then
                print *, "Number of parcels disagree!"
                stop
            endif
        end subroutine init_regular_positions

        subroutine init_refine(s_parcels, n_par, lam)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            double precision,                    intent(inout) :: lam
            double precision                                   :: a2

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call do_ellipse_split(s_parcels, n_par, parcel%lambda_max)
                a2 = get_eigenvalue(s_parcels%B(:, 1))
                lam = a2 / get_ab(s_parcels%area(1))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of bilinear
        ! interpolation and "apar"
        subroutine alloc_and_precompute(s_parcels, n_par)
            type(surface_parcel_container_type), intent(in) :: s_parcels
            integer,                             intent(in) :: n_par
            double precision, allocatable                   :: resi(:, :)
            double precision                                :: rsum
            integer                                         :: n, l

            allocate(apar(n_par))
            allocate(weights(ngp, n_par))
            allocate(is(ngp, n_par))
            allocate(js(ngp, n_par))
            allocate(resi(0:ny-1, 0:nx-1))

            ! Compute mean parcel density:
            resi = zero

            !$omp parallel do default(shared) private(n, l) reduction(+:resi)
            do n = 1, n_par
                ! get interpolation weights and mesh indices
                call bilinear(s_parcels%position(:, n), is(:, n), js(:, n), weights(:, n))

                do l = 1, ngp
                    resi(js(l, n), is(l, n)) = resi(js(l, n), is(l, n)) + weights(l, n)
                enddo
            enddo
            !$omp end parallel do

            ! Determine local inverse density of parcels (apar)
            !$omp parallel do default(shared) private(n, l, rsum)
            do n = 1, n_par
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
        subroutine init_from_grids(s_parcels, n_par, which, ncfname, tol)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            character(2),                        intent(in)    :: which
            character(*),                        intent(in)    :: ncfname
            double precision,                    intent(in)    :: tol
            double precision                                   :: buffer(0:nz, 0:ny-1, 0:nx-1)
            integer                                            :: ncid, k
            integer                                            :: n_steps, start(4), cnt(4)

            call alloc_and_precompute(s_parcels, n_par)

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            cnt  =  (/ nx, ny, nz+1, 1       /)
            start = (/ 1,  1,  1,    n_steps /)

            k = nz
            if (which == 'lo') then
                k = 0
            endif


            if (has_dataset(ncid, 'x_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'x_vorticity', buffer, start=start, cnt=cnt)
                call gen_parcel_scalar_attr(s_parcels, n_par, &
                                            buffer(k, :, :), tol, s_parcels%vorticity(1, :))
            endif

            if (has_dataset(ncid, 'y_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'y_vorticity', buffer, start=start, cnt=cnt)
                call gen_parcel_scalar_attr(s_parcels, n_par, &
                                            buffer(k, :, :), tol, s_parcels%vorticity(2, :))
            endif

            if (has_dataset(ncid, 'z_vorticity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'z_vorticity', buffer, start=start, cnt=cnt)
                call gen_parcel_scalar_attr(s_parcels, n_par, &
                                            buffer(k, :, :), tol, s_parcels%vorticity(3, :))
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'buoyancy', buffer, start=start, cnt=cnt)
                call gen_parcel_scalar_attr(s_parcels, n_par, &
                                            buffer(k, :, :), tol, s_parcels%buoyancy)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'humidity')) then
                buffer = zero
                call read_netcdf_dataset(ncid, 'humidity', buffer, start=start, cnt=cnt)
                call gen_parcel_scalar_attr(s_parcels, n_par, &
                                            buffer(k, :, :), tol, s_parcels%humidity)
            endif
#endif
            call close_netcdf_file(ncid)

            call dealloc

        end subroutine init_from_grids

        ! Generates the parcel attribute "par" from the field values provided
        ! in "field" (see Fontane & Dritschel, J. Comput. Phys. 2009, section 2.2)
        subroutine gen_parcel_scalar_attr(s_parcels, n_par, field, tol, par)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            double precision,                    intent(in)  :: field(0:ny-1, 0:nx-1)
            double precision,                    intent(in)  :: tol
            double precision,                    intent(out) :: par(:)
            double precision                                 :: resi(0:ny-1, 0:nx-1)
            double precision                                 :: rms, rtol, rerr
            double precision                                 :: rsum, fsum, avg_field
            integer                                          :: n, l

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print *, 'Generate parcel attribute'
                endif
#endif

            ! Compute mean field value:
            ! (divide by nx * ny since lower and upper edge weights are halved)
            avg_field = sum(field) / dble(nx * ny)

            resi = (field - avg_field) ** 2

            rms = dsqrt(sum(resi) / dble(nx * ny))


            if (rms == zero) then
                !$omp parallel default(shared)
                !$omp do private(n)
                do n = 1, n_par
                    ! assign mean value
                    par(n) = avg_field
                enddo
                !$omp end do
                !$omp end parallel
                return
            endif

            ! Maximum error permitted below in gridded residue:
            rtol = rms * tol

            ! Initialise (area-weighted) parcel attribute with a guess
            !$omp parallel do default(shared) private(n, l, fsum)
            do n = 1, n_par
                fsum = zero
                do l = 1, ngp
                    fsum = fsum + field(js(l, n), is(l, n)) * weights(l, n)
                enddo
                par(n) = apar(n) * fsum
            enddo
            !$omp end parallel do

            ! Iteratively compute a residual and update (area-weighted) attribute:
            rerr = one

            do while (rerr .gt. rtol)
                !Compute residual:
                resi = zero
                do n = 1, n_par
                    do l = 1, ngp
                        resi(js(l, n), is(l, n)) = resi(js(l, n), is(l, n)) &
                                                 + weights(l, n) * par(n)
                    enddo
                enddo

                resi = field - resi

                !Update (area-weighted) attribute:
                !$omp parallel do default(shared) private(n, rsum, l)
                do n = 1, n_par
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

            !Finally divide by parcel area to define attribute:
            ! (multiply with acell since algorithm is designed for area fractions)
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                par(n) = acell * par(n) / s_parcels%area(n)
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine gen_parcel_scalar_attr

end module surface_parcel_init
