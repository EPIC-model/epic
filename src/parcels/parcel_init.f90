! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parcel_interpl, only : trilinear, ngp
    use parameters, only : update_parameters,   &
                           dx, vcell, ncell,    &
                           extent, lower, nx, nz
    use h5_reader
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: init_handle

    double precision, allocatable :: weights(:, :), apar(:)
    integer, allocatable :: is(:, :), js(:, :)

    private :: weights, apar, is, js


    private :: init_regular_positions,      &
               init_refine,                 &
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


        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(h5fname, tol)
            character(*),     intent(in) :: h5fname
            double precision, intent(in) :: tol
            double precision             :: lam, ratio
            integer(hid_t)               :: h5handle

            call start_timer(init_handle)

            ! read domain dimensions
            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)
            call read_h5_box(h5handle, nx, nz, extent, lower)
            call close_h5_file(h5handle)

            ! update global parameters
            call update_parameters

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel%n_per_cell * ncell

            call init_regular_positions

            ! initialize the volume of each parcel
            parcels%volume(1:n_parcels) = vcell / dble(parcel%n_per_cell)

            if (parcel%is_elliptic) then
                deallocate(parcels%stretch)

                ratio = dx(1) / dx(2)

                ! aspect ratio: lam = a / b
                lam = max(dx(2) / dx(1), ratio)

                ! B11
                parcels%B(1:n_parcels, 1) = ratio * get_ab(parcels%volume(1:n_parcels))

                ! B12
                parcels%B(1:n_parcels, 2) = zero

                call init_refine(lam)

            else
                deallocate(parcels%B)
                parcels%stretch = zero
            endif

            parcels%vorticity(1:n_parcels) = zero
            parcels%buoyancy(1:n_parcels) = zero
            parcels%humidity(1:n_parcels) = zero

            call init_from_grids(h5fname, tol)

            call stop_timer(init_handle)

        end subroutine init_parcels


        subroutine init_regular_positions
            integer          :: i, ii, j, jj, k, n_per_dim
            double precision :: del(2)

            ! number of parcels per dimension
            n_per_dim = dsqrt(dble(parcel%n_per_cell))
            if (n_per_dim ** 2 .ne. parcel%n_per_cell) then
                print *, "Number of parcels per cell (", &
                         parcel%n_per_cell, ") not a square."
                stop
            endif

            del = dx / dble(two * n_per_dim)

            k = 1
            do j = 0, nz-1
                do i = 0, nx-1
                    do jj = 1, 2 * n_per_dim, 2
                        do ii = 1, 2 * n_per_dim, 2
                            parcels%position(k, 1) = lower(1) + dble(i) * dx(1) + del(1) * dble(ii)
                            parcels%position(k, 2) = lower(2) + dble(j) * dx(2) + del(2) * dble(jj)
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

            if (.not. parcel%is_elliptic) then
                return
            endif

            ! do refining by splitting
            do while (lam >= parcel%lambda)
                call split_ellipses(parcels, parcel%lambda, parcel%vmaxfraction)
                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine


        ! Precompute weights, indices of trilinear
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
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(n, :), is(n, :), js(n, :), weights(n, :))

                do l = 1, ngp
                    resi(js(n, l), is(n, l)) = resi(js(n, l), is(n, l)) + weights(n, l)
                enddo
            enddo

            !Double edge values at iz = 0 and nz:
            resi(0, :) = two * resi(0, :)
            resi(nz,:) = two * resi(nz,:)

            ! Determine local inverse density of parcels (apar)
            do n = 1, n_parcels
                rsum = zero
                do l = 1, ngp
                    rsum = rsum + resi(js(n, l), is(n, l)) * weights(n, l)
                enddo
                apar(n) = one / rsum
            enddo

        end subroutine alloc_and_precompute

        subroutine dealloc
            deallocate(apar)
            deallocate(weights)
            deallocate(is)
            deallocate(js)
        end subroutine dealloc


        subroutine init_from_grids(h5fname, tol)
            character(*),     intent(in)  :: h5fname
            double precision, intent(in)  :: tol
            double precision, allocatable :: buffer_2d(:, :)
            double precision              :: field_2d(-1:nz+1, 0:nx-1)
            integer(hid_t)                :: h5handle

            call alloc_and_precompute

            call open_h5_file(h5fname, H5F_ACC_RDONLY_F, h5handle)

            if (has_dataset(h5handle, 'vorticity')) then
                call read_h5_dataset_2d(h5handle, 'vorticity', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, tol, parcels%vorticity)
            endif


            if (has_dataset(h5handle, 'buoyancy')) then
                call read_h5_dataset_2d(h5handle, 'buoyancy', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, tol, parcels%buoyancy)
            endif

            call close_h5_file(h5handle)

            call dealloc

        end subroutine init_from_grids

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
                ! assign mean value
                par(1:n_parcels) = avg_field
                return
            endif

            ! Maximum error permitted below in gridded residue:
            rtol = rms * tol

            ! Initialise (volume-weighted) parcel attribute with a guess
            do n = 1, n_parcels
                fsum = zero
                do l = 1, ngp
                    fsum = fsum + field(js(n, l), is(n, l)) * weights(n, l)
                enddo
                par(n) = apar(n) * fsum
            enddo

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
                do n = 1, n_parcels
                    rsum = zero
                    do l = 1, ngp
                        rsum = rsum + resi(js(n, l), is(n, l)) * weights(n, l)
                    enddo
                    par(n) = par(n) + apar(n) * rsum
                enddo

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
            par(1:n_parcels) = vcell * par(1:n_parcels) / parcels%volume(1:n_parcels)

        end subroutine gen_parcel_scalar_attr

end module parcel_init
