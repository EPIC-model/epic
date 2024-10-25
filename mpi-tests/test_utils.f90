module test_utils
    use constants, only : zero, f12, f23, one, two, pi, twopi
    use parameters, only : lower, vmin, dx, nz, center
    use mpi_timer
    use parcel_container, only : resize_timer, n_parcels, parcels
    use parcel_split_mod, only : split_timer
    use parcel_merging, only : merge_timer
    use parcel_nearest, only : merge_nearest_timer          &
                             , merge_tree_resolve_timer     &
                             , nearest_allreduce_timer      &
                             , nearest_barrier_timer        &
                             , nearest_rma_timer
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
    use parcel_interpl, only : grid2par_timer, par2grid_timer, halo_swap_timer
    use parcel_init, only : init_timer
    use ls_rk, only : rk_timer
    use mpi_environment
    use mpi_layout
    use options, only : parcel
    use mpi_utils, only : mpi_exit_on_error
    use bndry_fluxes, only : bndry_flux_timer
    use parcel_damping, only : damping_timer
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
            call register_timer('MPI allreduce timer (in tree resolve)', nearest_allreduce_timer)
            call register_timer('MPI barrier timer (in tree resolve)', nearest_barrier_timer)
            call register_timer('MPI RMA timer (in tree resolve)', nearest_rma_timer)
            call register_timer('p2g/v2g halo (non-excl.)', halo_swap_timer)
            call register_timer('boundary fluxes', bndry_flux_timer)
            call register_timer('damping', damping_timer)
        end subroutine register_all_timers

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine init_rng

            call random_seed(size = sk)
            allocate(seed(1:sk))
            call random_seed(get=seed)

        end subroutine init_rng

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Box-Muller transform
        ! (5 April 2024, https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
        subroutine random_normal(u1, u2)
            double precision, intent(inout) :: u1, u2
            double precision                :: z

            z = sqrt(-two * log(u1))
            u1 = z * cos(twopi * u2)
            u2 = z * sin(twopi * u2)

        end subroutine random_normal

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ! pick point uniformly on a sphere:
        ! 5 April 2024, https://stats.stackexchange.com/a/7984
        subroutine random_angles(theta, phi)
            double precision, intent(out) :: theta, phi
            double precision              :: u(4)

            ! get 4 uniform numbers in [0, 1[:
            call random_number(u)

            ! transform to standard normal:
            call random_normal(u(1), u(2))
            call random_normal(u(3), u(4))

            ! normalise (note: we do not need u(4) later on):
            u(4) = u(1) ** 2 + u(2) ** 2 + u(3) ** 2
            u(1:3) = u(1:3) / u(4)

            ! azimuthal angle, [0, 2pi[
            theta = datan2(u(2), u(1))

            ! polar angle, [0, pi[
            u(3) = max(-1.0d0, min(u(3), 1.0d0))
            phi = dacos(u(3))

        end subroutine random_angles

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine setup_parcels(xlen, ylen, zlen, l_shuffle, l_variable_nppc)
            double precision, intent(in) :: xlen, ylen, zlen
            logical,          intent(in) :: l_shuffle, l_variable_nppc
            double precision             :: rn(10), lam, lam2, abc, a2, b2, c2, theta, phi
            double precision             :: st, ct, sp, cp, corner(3), xhw, yhw, zhw, x, y, z
            double precision             :: xlo, xhi, ylo, yhi, zlo, zhi
            integer                      :: ix, iy, iz, m, l, npp, n_per_cell

            npp = parcel%n_per_cell
            n_per_cell = max(10, parcel%n_per_cell - 10)

            xhw = f12 * xlen
            xlo = center(1) - xhw
            xhi = xlo + xlen

            yhw = f12 * ylen
            ylo = center(2) - yhw
            yhi = ylo + ylen

            zhw = f12 * zlen
            zlo = center(3) - zhw
            zhi = zlo + zlen

            l = 1
            do iz = 0, nz-1
                do iy = box%lo(2), box%hi(2)
                   do ix = box%lo(1), box%hi(1)
                        if (l_variable_nppc) then
                            call random_number(lam)
                            npp = int(lam * n_per_cell) + 10    ! ensure at least 10 parcels per cell
                        endif
                        corner = lower + dble((/ix, iy, iz/)) * dx
                        do m = 1, npp
                            ! rn between 0 and 1
                            call random_number(rn)

                            x = corner(1) + dx(1) * rn(1)
                            y = corner(2) + dx(2) * rn(2)
                            z = corner(3) + dx(3) * rn(3)

                            parcels%position(1, l) = x
                            parcels%position(2, l) = y
                            parcels%position(3, l) = z

                            ! vorticity between -10 and 10: y = 20 * x - 10
                            parcels%vorticity(1, l) = 20.0d0 * rn(4) - 10.d0
                            parcels%vorticity(2, l) = 20.0d0 * rn(5) - 10.d0
                            parcels%vorticity(3, l) = 20.0d0 * rn(6) - 10.d0

                            ! buoyancy between -1 and 1: y = 2 * x - 1
                            parcels%buoyancy(l) = 2.0d0 * rn(7) - 1.d0

                            if ((x >= xlo) .and. (x <= xhi)  .and. &
                                (y >= ylo) .and. (y <= yhi)  .and. &
                                (z >= zlo) .and. (z <= zhi)) then

                                ! volume between 0.5 * vmin and 1.5 * vmin
                                parcels%volume(l) = vmin * rn(8) + f12 * vmin
                            else
                                ! volume between vmin and 2 * vmin
                                parcels%volume(l) = vmin * rn(8) + vmin
                            endif


                            ! lam = a / c in [1, 4]
                            lam = 3.d0 * rn(9) + 1.0d0

                            ! lam2 = a / b
                            lam2 = 3.d0 * rn(10) + 1.0d0

                            abc = 0.75d0 * parcels%volume(l) / pi

                            a2 = (abc * lam * lam2)  ** f23
                            b2 = a2 / lam2 ** 2
                            c2 = a2 / lam ** 2

                            ! get random angles:
                            call random_angles(theta, phi)

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
             n_parcels = l - 1

            if (l_shuffle) then
                call shuffleall
            endif

        end subroutine setup_parcels

        ! performs a Fisher-Yates shuffle (aka Knuth shuffle)
        subroutine shuffleall
            integer          :: shuffle_index, rand_target
            double precision :: tmp_var
            double precision :: tmp_vec(3), tmp_B(5)
            double precision :: random_out

            do shuffle_index = n_parcels, 2, -1
               call random_number(random_out)
               rand_target = int(random_out * shuffle_index) + 1

               tmp_vec = parcels%position(:, rand_target)
               parcels%position(:, rand_target) = parcels%position(:, shuffle_index)
               parcels%position(:, shuffle_index) = tmp_vec

               tmp_vec = parcels%vorticity(:, rand_target)
               parcels%vorticity(:, rand_target) = parcels%vorticity(:, shuffle_index)
               parcels%vorticity(:, shuffle_index) = tmp_vec

               tmp_var = parcels%buoyancy(rand_target)
               parcels%buoyancy(rand_target) = parcels%buoyancy(shuffle_index)
               parcels%buoyancy(shuffle_index) = tmp_var

               tmp_var = parcels%volume(rand_target)
               parcels%volume(rand_target) = parcels%volume(shuffle_index)
               parcels%volume(shuffle_index) = tmp_var

               tmp_B = parcels%B(:, rand_target)
               parcels%B(:, rand_target) = parcels%B(:, shuffle_index)
               parcels%B(:, shuffle_index) = tmp_B
            end do
        end subroutine shuffleall

end module test_utils
