! =============================================================================
!               This module initializes parcels for the static
!               Taylor-Green vortex test case.
! =============================================================================
module parcel_tg_init
    use options, only : parcel
    use constants, only : zero, two, one, f12, f34, pi, twopi, max_num_parcels
    use parcel_container, only : parcels, n_parcels
    use parcel_ellipse, only : get_ab, get_B22, get_eigenvalue
    use parcel_split, only : split_ellipses
    use parameters, only : update_parameters,   &
                           dx, vcell, ncell,    &
                           extent, lower, nx, nz
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none

    integer :: init_timer

    private :: init_refine

    contains

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_tg_parcels
            double precision             :: lam, ratio, dr, dr2, rpos2, x0, z0, xp, zp
            integer                      :: n

            call start_timer(init_timer)

            lower = zero
            extent = one
            nx = 100
            nz = 100
            ! disk centre (x0, z0) and radius (dr)
            dr = 0.15d0
            x0 = f12
            z0 = f34

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
                parcels%B(1, n) = ratio * get_ab(parcels%volume(n))

                ! B12
                parcels%B(2, n) = zero
            enddo
            !$omp end do
            !$omp end parallel

            call init_refine(lam)


            dr2 = dr ** 2
            !$omp parallel default(shared)
            !$omp do private(n, rpos2, xp, zp)
            do n = 1, n_parcels
                parcels%buoyancy(n) = zero

                ! humidity is tracer
                xp = parcels%position(1, n)
                zp = parcels%position(2, n)
                rpos2 = (xp - x0) ** 2 + (zp - z0) ** 2

                parcels%humidity(n) = one
                if (rpos2 <= dr2) then
                    parcels%humidity(n) = two
                endif

                ! TG vorticity
                parcels%vorticity(n) = twopi * dsin(pi * xp) * dsin(pi * zp)
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(init_timer)

        end subroutine init_tg_parcels


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
                call split_ellipses(parcels, parcel%lambda_max)
                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine init_refine


end module parcel_tg_init
