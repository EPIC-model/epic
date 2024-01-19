! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parcel_interpl, only : bilinear, ngp
    use fields, only : vortg, tbuoyg
    use parameters, only : dx, vcell, ncell,        &
                           extent, lower, nx, nz,   &
                           max_num_parcels
    use netcdf_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    use surface_parcel_init, only : init_surface_parcels
    implicit none

    integer :: init_timer

    private :: init_refine, init_parcels_from_grids

    contains

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_parcels(fname)
            character(*),     intent(in) :: fname
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

            call init_parcels_from_grids(fname)

            call init_surface_parcels(fname)

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
            double precision                :: B22, a2

            ! do refining by splitting
            do while (lam >= parcel%lambda_max)
                call split_ellipses
                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_parcels_from_grids(ncfname)
            character(*),     intent(in)  :: ncfname
            integer                        :: n, l
            integer                       :: ncid
            integer                       :: n_steps, start(3), cnt(3)
            integer :: ii(4), jj(4)
            double precision :: ww(4)

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            cnt  =  (/ nx, nz+1, 1       /)
            start = (/ 1,  1,    n_steps /)

            vortg = zero
            tbuoyg = zero

            if (has_dataset(ncid, 'vorticity')) then
                call read_netcdf_dataset(ncid, 'vorticity', vortg(0:nz, :), start=start, cnt=cnt)
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                call read_netcdf_dataset(ncid, 'buoyancy', tbuoyg(0:nz, :), start=start, cnt=cnt)
            endif

            !$omp parallel default(shared)
            !$omp do private(n, l, ii, jj, ww)
            do n = 1, n_parcels

                ! get interpolation weights and mesh indices
                call bilinear(parcels%position(:, n), ii, jj, ww)

                ! loop over grid points which are part of the interpolation
                do l = 1, 4
                    parcels%vorticity(n) = parcels%vorticity(n) + ww(l) * vortg(jj(l), ii(l))
                parcels%buoyancy(n) = parcels%buoyancy(n) + ww(l) * tbuoyg(jj(l), ii(l))
                enddo
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(init_timer)

        end subroutine init_parcels_from_grids

end module parcel_init
