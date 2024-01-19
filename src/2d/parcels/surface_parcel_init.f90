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
    use parameters, only : dx, lcell, vcell,        &
                           extent, lower, nx, nz,   &
                           max_num_surf_parcels
    use netcdf_reader
    use fields, only : tbuoyg, vortg
    use timer, only : start_timer, stop_timer
    use omp_lib
    implicit none


    private :: init_from_grids,             &
               init_surface_parcels_

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
            integer                                            :: n, n_per_dim
            double precision                                   :: length

            n_per_dim = int(dsqrt(dble(parcel%n_per_cell)))

            ! we use "n_per_cell" parcels per grid cell
            n_par = n_per_dim * nx

            if (n_par > max_num_surf_parcels) then
                print *, "Number of parcels exceeds limit of", &
                          max_num_surf_parcels, ". Exiting."
                stop
            endif

            !------------------------------------------------------------------
            ! set up parcel positions and 'right' array:
            length = lcell / dble(n_per_dim)

            spar%position(1) = lower(1)
            spar%right(1) = 2
            do n = 2, n_par
                spar%position(n) = spar%position(n-1) + length
                spar%right(n) = mod(n, n_par) + 1
            enddo
            spar%volume(1:n_par) = vcell / dble(2.0d0)


            !------------------------------------------------------------------
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
