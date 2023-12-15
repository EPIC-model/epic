! =============================================================================
!               This module initializes surface parcel default values.
! =============================================================================
module surface_parcel_init
    use options, only : parcel, output, verbose, field_tol
    use constants, only : zero, two, one, f12
    use surface_parcel_container, only : n_top_parcels, top_parcels &
                                       , n_bot_parcels, bot_parcels &
                                       , surface_parcel_container_type
    use surface_parcel_interpl, only : linear, ngp
    use parameters, only : dx, lcell,               &
                           extent, lower, nx, nz,   &
                           max_num_surf_parcels
    use netcdf_reader
    use fields, only : tbuoyg, vortg
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none


    private :: init_from_grids,             &
               init_surface_parcels_,       &
               init_regular_positions_

    contains

        subroutine init_surface_parcels(fname)
            character(*), intent(in) :: fname

            call init_surface_parcels_(0,  n_bot_parcels, bot_parcels, fname)
            call init_surface_parcels_(nz, n_top_parcels, top_parcels, fname)

        end subroutine init_surface_parcels

        ! Set default values for parcel attributes
        ! Attention: This subroutine assumes that the parcel
        !            container is already allocated!
        subroutine init_surface_parcels_(iz, n_par, spar, fname)
            integer,                             intent(in)    :: iz
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            character(*),                        intent(in)    :: fname
            integer                                            :: n

            ! we use "n_per_cell" parcels per grid cell
            n_par = parcel%n_per_cell * nx

            if (n_par > max_num_surf_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_surf_parcels, ". Exiting."
                stop
            endif

            call init_regular_positions_(n_par, spar)

            ! initialize the length of each parcel
            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                spar%length(n) = lcell / dble(parcel%n_per_cell)
            enddo
            !$omp end do
            !$omp end parallel

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_par
                spar%vorticity(n) = zero
                spar%buoyancy(n) = zero
#ifndef ENABLE_DRY_MODE
                spar%humidity(n) = zero
#endif
            enddo
            !$omp end do
            !$omp end parallel

            call init_from_grids(iz, n_par, spar, fname)

        end subroutine init_surface_parcels_


        ! Position parcels regularly in the domain.
        subroutine init_regular_positions_(n_par, spar)
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer                                            :: ix, i, k
            double precision                                   :: im, corner

            im = one / dble(parcel%n_per_cell)

            k = 1
            do ix = 0, nx-1
                corner = lower(1) + dble(ix) * dx(1)
                do i = 1, parcel%n_per_cell
                    spar%position(k) = corner + dx(1) * (dble(i) - f12) * im
                    k = k + 1
                enddo
            enddo

            if (.not. n_par == k - 1) then
                print *, "Number of surface parcels disagree!"
                stop
            endif
        end subroutine init_regular_positions_

        ! Initialise parcel attributes from gridded quantities.
        ! Attention: This subroutine currently only supports
        !            vorticity and buoyancy fields.
        subroutine init_from_grids(iz, n_par, spar, ncfname)
            integer,                             intent(in)    :: iz
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            character(*),                        intent(in)    :: ncfname
            integer                                            :: ncid
            integer                                            :: n_steps, start(3), cnt(3)
            double precision                                   :: weights(2)
            integer                                            :: is(2), l, n

            call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)

            call get_num_steps(ncid, n_steps)

            cnt  =  (/ nx, nz+1, 1       /)
            start = (/ 1,  1,    n_steps /)

            vortg  = zero
            tbuoyg = zero

            if (has_dataset(ncid, 'vorticity')) then
                call read_netcdf_dataset(ncid, 'vorticity', vortg(0:nz, :), start=start, cnt=cnt)
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                call read_netcdf_dataset(ncid, 'buoyancy', tbuoyg(0:nz, :), start=start, cnt=cnt)
            endif

            call close_netcdf_file(ncid)

            !$omp parallel default(shared)
            !$omp do private(n, l, is, weights)
            do n = 1, n_par

                ! get interpolation weights and mesh indices
                call linear(spar%position(n), is, weights)

                ! loop over grid points which are part of the interpolation
                do l = 1, 2
                    spar%vorticity(n) = spar%vorticity(n) + weights(l) * vortg(iz, is(l))
                    spar%buoyancy(n)  = spar%buoyancy(n) + weights(l) * tbuoyg(iz, is(l))
                enddo
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine init_from_grids

end module surface_parcel_init
