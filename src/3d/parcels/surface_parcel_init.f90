! =============================================================================
!             This module initializes surface parcel default values.
! =============================================================================
module surface_parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, one, f12
    use surface_parcel_container, only : surface_parcel_container_type      &
                                       , lo_surf_parcels, n_lo_surf_parcels &
                                       , up_surf_parcels, n_up_surf_parcels &
                                       , surface_parcel_alloc
    use parcel_ellipse, only : get_ab, get_eigenvalue
    use surface_parcel_split_mod, only : do_ellipse_split
    use surface_parcel_interpl, only : bilinear, ngp
    use parameters, only : dx, acell,                   &
                           extent, lower, nx, ny, nz,   &
                           max_num_surf_parcels
    use netcdf_reader
    use timer, only : start_timer, stop_timer
    use omp_lib
    use fields, only : vortg, tbuoyg
#ifndef ENABLE_DRY_MODE
    use fields, only : humg
#endif
    implicit none

    integer :: surf_init_timer

    double precision :: weights(ngp)
    integer :: is(ngp), js(ngp)

    private :: weights, is, js


    private :: init_refine,                 &
               parcel_default,              &
               init_from_grids

    contains

        subroutine surface_parcel_default
            call parcel_default(lo_surf_parcels, n_lo_surf_parcels)
            call parcel_default(up_surf_parcels, n_up_surf_parcels)
        end subroutine surface_parcel_default

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine parcel_default(s_parcels, n_par)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
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

            call surface_parcel_alloc(s_parcels, max_num_surf_parcels)


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

            call stop_timer(surf_init_timer)

        end subroutine parcel_default


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
                call do_ellipse_split(s_parcels, n_par)
                a2 = get_eigenvalue(s_parcels%B(:, 1))
                lam = a2 / get_ab(s_parcels%area(1))
            end do
        end subroutine init_refine


        subroutine init_surface_from_grids
            call init_from_grids(lo_surf_parcels, n_lo_surf_parcels, 'lo')
            call init_from_grids(up_surf_parcels, n_up_surf_parcels, 'up')
        end subroutine init_surface_from_grids

        ! Initialise parcel attributes from gridded quantities.
        subroutine init_from_grids(s_parcels, n_par, which)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            character(2),                        intent(in)    :: which
            integer                                            :: n, k, l

            k = nz
            if (which == 'lo') then
                k = 0
            endif

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, weights)
            do n = 1, n_par

                ! get interpolation weights and mesh indices
                call bilinear(s_parcels%position(:, n), is, js, weights)

                ! loop over grid points which are part of the interpolation
                do l = 1, ngp
                    s_parcels%vorticity(:, n) = s_parcels%vorticity(:, n) &
                                              + weights(l) * vortg(k, js(l), is(l), :)
                    s_parcels%buoyancy(n) = s_parcels%buoyancy(n) &
                                          + weights(l) * tbuoyg(k, js(l), is(l))
#ifndef ENABLE_DRY_MODE
                    s_parcels%humidity(n) = s_parcels%humidity(n) &
                                          + weights(l) * humg(k, js(l), is(l))
#endif

                enddo
            enddo
            !$omp end do
            !$omp end parallel
        end subroutine init_from_grids

end module surface_parcel_init
