module test_utils
    use constants, only : pi, zero, one, f12, f23, twopi
    use parameters, only : lower, vmin, dx, nz
    use mpi_timer
    use parcel_container, only : resize_timer, n_parcels, parcels
    use parcel_split_mod, only : split_timer
    use parcel_merging, only : merge_timer
    use parcel_nearest, only : merge_nearest_timer, merge_tree_resolve_timer
    use parcel_correction, only : lapl_corr_timer,        &
                                  grad_corr_timer,        &
                                  vort_corr_timer
    use parcel_diagnostics, only : parcel_stats_timer
    use parcel_netcdf, only : parcel_io_timer
    use parcel_diagnostics_netcdf, only : parcel_stats_io_timer
    use field_netcdf, only : field_io_timer
    use field_diagnostics, only : field_stats_timer
    use field_diagnostics_netcdf, only : field_stats_io_timer
    use inversion_mod, only : vor2vel_timer, vtend_timer
    use parcel_interpl, only : grid2par_timer, par2grid_timer
    use parcel_init, only : init_timer
    use ls_rk, only : rk_timer
    use mpi_environment
    use mpi_layout
    use options, only : parcel
    use mpi_utils, only : mpi_exit_on_error
    implicit none

    integer :: epic_timer
    integer, allocatable :: seed(:)
    integer :: sk


    contains

        subroutine register_all_timers
            call register_timer('epic', epic_timer)
            call register_timer('parcel container resize', resize_timer)
            call register_timer('par2grid', par2grid_timer)
            call register_timer('grid2par', grid2par_timer)
            call register_timer('parcel split', split_timer)
            call register_timer('parcel merge', merge_timer)
            call register_timer('laplace correction', lapl_corr_timer)
            call register_timer('gradient correction', grad_corr_timer)
            call register_timer('net vorticity correction', vort_corr_timer)
            call register_timer('parcel initialisation', init_timer)
            call register_timer('parcel diagnostics', parcel_stats_timer)
            call register_timer('parcel I/O', parcel_io_timer)
            call register_timer('parcel diagnostics I/O', parcel_stats_io_timer)
            call register_timer('field I/O', field_io_timer)
            call register_timer('field diagnostics', field_stats_timer)
            call register_timer('field diagnostics I/O', field_stats_io_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('parcel push', rk_timer)
            call register_timer('merge nearest', merge_nearest_timer)
            call register_timer('merge tree resolve', merge_tree_resolve_timer)
        end subroutine register_all_timers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_rng
            call random_seed(size=sk)
            allocate(seed(1:sk))
            seed(:) = world%rank
            call random_seed(put=seed)
        end subroutine init_rng

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine setup_parcels
            double precision :: rn(12), lam, lam2, abc, a2, b2, c2, theta, phi
            double precision :: st, ct, sp, cp, corner(3)
            integer          :: ix, iy, iz, m, l

            n_parcels = parcel%n_per_cell * box%ncell

            l = 1
            do iz = 0, nz-1
                do iy = box%lo(2), box%hi(2)
                    do ix = box%lo(1), box%hi(1)
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do m = 1, parcel%n_per_cell
                            ! rn between 0 and 1
                            call random_number(rn)

                            parcels%position(1, l) = corner(1) + dx(1) * rn(1)
                            parcels%position(2, l) = corner(2) + dx(2) * rn(2)
                            parcels%position(3, l) = corner(3) + dx(3) * rn(3)

                            ! vorticity between -10 and 10: y = 20 * x - 10
                            parcels%vorticity(1, l) = 20.0d0 * rn(4) - 10.d0
                            parcels%vorticity(2, l) = 20.0d0 * rn(5) - 10.d0
                            parcels%vorticity(3, l) = 20.0d0 * rn(6) - 10.d0

                            ! buoyancy between -1 and 1: y = 2 * x - 1
                            parcels%buoyancy(l) = 2.0d0 * rn(7) - 1.d0

                            ! volume between 0.5 * vmin and 1.5 * vmin
                            parcels%volume(l) = vmin * rn(8) + f12 * vmin

                            ! lam = a / c in [1, 4]
                            lam = 3.d0 * rn(9) + 1.0d0

                            ! lam2 = a / b
                            lam2 = 3.d0 * rn(10) + 1.0d0

                            abc = 0.75d0 * parcels%volume(l) / pi

                            a2 = (abc * lam * lam2)  ** f23
                            b2 = a2 / lam2 ** 2
                            c2 = a2 / lam ** 2

                            ! theta in [0, 2pi[
                            theta = twopi * rn(11)

                            ! phi in [0, pi[
                            phi = pi * rn(12)

                            st = dsin(theta)
                            ct = dcos(theta)
                            sp = dsin(phi)
                            cp = dcos(phi)

                            parcels%B(1, l) = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
                            parcels%B(2, l) = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
                            parcels%B(3, l) = (a2 - c2) * ct * sp * cp
                            parcels%B(4, l) = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
                            parcels%B(5, l) = (a2 - c2) * st * sp * cp

                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo

            if (.not. n_parcels == l - 1) then
                call mpi_exit_on_error("Number of parcels disagree!")
            endif

        end subroutine setup_parcels

end module test_utils
