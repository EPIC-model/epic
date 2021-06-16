! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output
    use constants, only : zero, two, one, f12
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parcel_interpl, only : trilinear, ngp
    use parameters, only : update_parameters,       &
                           write_h5_parameters,     &
                           dx, vcell, ncell         &
                           extent, lower, nx, nz
    use reader
    implicit none

    private :: init_random_positions,       &
               init_regular_positions,      &
               init_refine,                 &
               init_from_grids,             &
               fill_field_from_buffer_2d

    contains


        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(filename)
            character(*), intent(in) :: filename
            double precision         :: lam, ratio

            ! read domain dimensions
            call open_h5_file(filename)
            call read_h5_box(nx, nz, extent, lower)
            call close_h5_file

            ! update global parameters
            call update_parameters

            call write_h5_parameters(trim(output%h5fname))

            ! set the number of parcels (see parcels.f90)
            ! we use "n_per_cell" parcels per grid cell
            n_parcels = parcel%n_per_cell * ncell

            if (parcel%is_random) then
                call init_random_positions
            else
                call init_regular_positions
            endif

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

            parcels%velocity(1:n_parcels, :) = zero
            parcels%vorticity(1:n_parcels) = zero
            parcels%buoyancy(1:n_parcels) = zero
            parcels%humidity(1:n_parcels) = zero

            call init_from_grids(filename)

        end subroutine init_parcels


        subroutine init_random_positions
            double precision :: val
            integer          :: n, k
            integer, allocatable :: seed(:)

            call random_seed(size=k)
            allocate(seed(1:k))
            seed(:) = parcel%seed
            call random_seed(put=seed)

            do n = 1, n_parcels
                call random_number(val)
                parcels%position(n, 1)= lower(1) + val * extent(1)
                call random_number(val)
                parcels%position(n, 2) = lower(2) + val * extent(2)
            enddo

            deallocate(seed)
        end subroutine init_random_positions


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


        subroutine init_from_grids(filename)
            character(*), intent(in)      :: filename
            double precision, allocatable :: buffer_2d(:, :)
            double precision              :: field_2d(-1:nz+1, 0:nx-1)

            call open_h5_file(filename)

            if (has_dataset('vorticity')) then
                call read_h5_dataset_2d('vorticity', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, 1.0d-9, parcels%vorticity)
            endif


            if (has_dataset('buoyancy')) then
                call read_h5_dataset_2d('buoyancy', buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, field_2d)
                deallocate(buffer_2d)
                call gen_parcel_scalar_attr(field_2d, 1.0d-9, parcels%buoyancy)
            endif

            call close_h5_file

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
            double precision :: apar(n_parcels)
            double precision :: rtol, rerr, rsum, fsum, fmean
            integer          :: is(ngp), js(ngp), n, l
            double precision :: weights(ngp)

            ! Compute mean field value:
            ! (divide by ncell since lower and upper edge weights are halved)
            fmean = (f12 * sum(field(0, :) + field(nz, :)) &
                         + sum(field(1:nz-1,:))) / dble(ncell)

            ! Maximum error permitted below in gridded residue:
            rtol = dabs(fmean) * tol

            if (rtol == zero) then
                ! keep initialisation to zero
                return
            endif

            ! Compute mean parcel density:
            resi = zero
            do n = 1, n_parcels

                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(n, :), is, js, weights)

                do l = 1, ngp
                    resi(js(l), is(l)) = resi(js(l), is(l)) + weights(l)
                enddo
            enddo

            !Double edge values at iz = 0 and nz:
            resi(0, :) = two * resi(0, :)
            resi(nz,:) = two * resi(nz,:)


            ! Determine local inverse density of parcels (apar)
            ! and initialise (volume-weighted) parcel attribute with a guess
            do n = 1, n_parcels
                ! get interpolation weights and mesh indices
                call trilinear(parcels%position(n, :), is, js, weights)

                rsum = zero
                fsum = zero
                do l = 1, ngp
                    rsum = rsum + resi(js(l), is(l)) * weights(l)
                    fsum = fsum + field(js(l), is(l)) * weights(l)
                enddo
                apar(n) = one / rsum
                par(n) = apar(n) * fsum
            enddo

            ! Iteratively compute a residual and update (volume-weighted) attribute:
            rerr = one

            do while (rerr .gt. rtol)
                !Compute residual:
                resi(0:nz, :) = zero
                do n = 1, n_parcels
                    ! get interpolation weights and mesh indices
                    call trilinear(parcels%position(n, :), is, js, weights)

                    do l = 1, ngp
                        resi(js(l), is(l)) = resi(js(l), is(l)) + weights(l) * par(n)
                    enddo
                enddo

                resi(0, :)    = two * resi(0, :)
                resi(nz, :)   = two * resi(nz, :)
                resi(0:nz, :) = field(0:nz, :) - resi(0:nz, :)

                !Update (volume-weighted) attribute:
                do n = 1, n_parcels
                    ! get interpolation weights and mesh indices
                    call trilinear(parcels%position(n, :), is, js, weights)

                    rsum = zero
                    do l = 1, ngp
                        rsum = rsum + resi(js(l), is(l)) * weights(l)
                    enddo
                    par(n) = par(n) + apar(n) * rsum
                enddo

                !Compute maximum error:
                rerr = maxval(dabs(resi))

                write(*,*) ' Max abs error = ',rerr
            enddo

            !Finally divide by parcel volume to define attribute:
            par(1:n_parcels) = par(1:n_parcels) / parcels%volume(1:n_parcels)

        end subroutine gen_parcel_scalar_attr

end module parcel_init
