!==============================================================================
! This module is used to apply boundary surface fluxes of buoyancy and
! humidity to each parcel in the lowest grid cell.
!==============================================================================
module bndry_fluxes
    use constants, only : zero, two, f13
    use parameters, only : dx, lower, dxi, nx, ny
    use parcel_interpl, only : bilinear
    use parcel_container, only : n_parcels, parcels
    use netcdf_reader
    use omp_lib
    use options, only : time
    use mpi_utils, only : mpi_stop
    use field_mpi, only : field_mpi_alloc                   &
                        , field_mpi_dealloc                 &
                        , field_buffer_to_halo              &
                        , field_interior_to_buffer          &
                        , interior_to_halo_communication
    use mpi_layout, only : box
    use mpi_timer, only : start_timer, stop_timer
    implicit none

    private

    logical, protected :: l_enable_flux
    logical            :: l_bndry_flux_allocated = .false.
    integer            :: bndry_flux_timer

    ! Spatial form of the buoyancy and humidity fluxes through lower surface:
    double precision, dimension(:, :), allocatable :: binc
#ifndef ENABLE_DRY_MODE
    double precision, dimension(:, :), allocatable :: hinc
#endif

    ! Denotes height below which surface fluxes are applied:
    double precision :: zdepth

    public :: l_enable_flux             &
            , apply_bndry_fluxes        &
            , read_bndry_fluxes         &
            , bndry_fluxes_allocate     &
            , bndry_fluxes_deallocate   &
            , bndry_flux_timer          &
            , bndry_fluxes_time_step

    contains

        subroutine bndry_fluxes_allocate
            if (l_bndry_flux_allocated) then
                return
            endif

            l_bndry_flux_allocated = .true.

            allocate(binc(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
#ifndef ENABLE_DRY_MODE
            allocate(hinc(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
#endif
        end subroutine bndry_fluxes_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_deallocate
            if (.not. l_bndry_flux_allocated) then
                return
            endif

            call start_timer(bndry_flux_timer)

            l_bndry_flux_allocated = .false.

            deallocate(binc)
#ifndef ENABLE_DRY_MODE
            deallocate(hinc)
#endif
            call stop_timer(bndry_flux_timer)

        end subroutine bndry_fluxes_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_default

            call start_timer(bndry_flux_timer)

            call bndry_fluxes_allocate

            binc = zero
#ifndef ENABLE_DRY_MODE
            hinc = zero
#endif
            call stop_timer(bndry_flux_timer)

        end subroutine bndry_fluxes_default

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Read in buoyancy and humidity fluxes.
        ! @pre The buoyancy boundary fluxes are assumed in units of m**2/s**3.
        subroutine read_bndry_fluxes(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, start(3), cnt(3)
            integer                  :: lo(3), hi(3)

            call start_timer(bndry_flux_timer)

            l_enable_flux = (fname /= '')

            if (.not. l_enable_flux) then
                return
            endif

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            lo = box%lo
            hi = box%hi

            ! need to add 1 since start must begin with index 1
            start(1:2) = lo(1:2) + 1
            start(3) = 1
            cnt(1:2) = hi(1:2) - lo(1:2) + 1
            cnt(3) = 1

            call bndry_fluxes_default

            if (has_dataset(ncid, 'bflux')) then
                call read_netcdf_dataset(ncid,                  &
                                         'bflux',               &
                                         binc(lo(2):hi(2),      &
                                              lo(1):hi(1)),     &
                                         start,                 &
                                         cnt)
            else
                call mpi_stop("No buoyancy flux field 'bflux' found in file.")
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'hflux')) then
                call read_netcdf_dataset(ncid,                  &
                                         'hflux',               &
                                         hinc(lo(2):hi(2),      &
                                              lo(1):hi(1)),     &
                                         start,                 &
                                         cnt)
            else
                call mpi_stop("No humidity flux field 'hflux' found in file.")
            endif
#endif

            call close_netcdf_file(ncid)


            call bndry_fluxes_fill_halo

            ! There is linear decrease of the buoyancy flux "force", where
            ! for z >= zdepth, the bouyancy fluy is zero.
            zdepth = lower(3) + dx(3)

            !------------------------------------------------------------------
            ! change the buoyancy flux units from m**2/s**3 to m/s**3;
            ! the time will be fixed when applying the flux
            ! Multiply with 2 to ensure that the flux given to the parcels
            ! entering the lowest grid layer matches the expected integrated flux:

            !$omp parallel
            !$omp workshare
            binc = two * dxi(3) * binc
#ifndef ENABLE_DRY_MODE
            hinc = two * dxi(3) * hinc
#endif
            !$omp end workshare
            !$omp end parallel

            call stop_timer(bndry_flux_timer)

        end subroutine read_bndry_fluxes

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Correct the time step estimate including the buoyancy flux.
        ! Uses g*dz/bf (where [bf] = m**2/s**3) to get a time measure.
        subroutine bndry_fluxes_time_step(dt)
            double precision, intent(inout) :: dt
            double precision                :: abs_max

            if (.not. l_enable_flux) then
                return
            endif

            ! local maximum of absolute value (units: m/s**3)
            abs_max = maxval(dabs(binc(box%lo(2):box%hi(2), box%lo(1):box%hi(1))))

            ! get global abs_max
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               abs_max,                 &
                               1,                       &
                               MPI_DOUBLE_PRECISION,    &
                               MPI_MAX,                 &
                               world%comm,              &
                               world%err)

            abs_max = (dx(3) / abs_max) ** f13

            dt = min(dt, time%alpha * abs_max)

        end subroutine bndry_fluxes_time_step

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Add fluxes of buoyancy and humidity to parcels in lowest grid layer:
        subroutine apply_bndry_fluxes(dt)
            double precision, intent(in) :: dt
            double precision             :: xy(2), z, fac, weights(0:1, 0:1)
            integer                      :: n, is, js

            if (.not. l_enable_flux) then
                return
            endif

            call start_timer(bndry_flux_timer)

            !$omp parallel default(shared)
            !$omp do private(n, is, js, weights, xy, z, fac)
            do n = 1, n_parcels
                z = parcels%position(3, n)
                if (z < zdepth) then

                    xy = parcels%position(1:2, n)

                    call bilinear(xy, is, js, weights)

                    ! "fac" has unit of time
                    ! Note: fac = 0 if z = zdepth
                    !       fac = dt if z = lower(3)
                    fac = (zdepth - z) * dxi(3) * dt

                    ! The multiplication by dt is necessary to provide the amount of b or h
                    ! entering through the bottom surface over a time interval of dt.
                    parcels%buoyancy(n) = parcels%buoyancy(n) &
                                        + fac * sum(weights * binc(js:js+1, is:is+1))
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n) = parcels%humidity(n) &
                                        + fac * sum(weights * hinc(js:js+1, is:is+1))
#endif
                endif
            enddo
            !$omp end do
            !$omp end parallel

            call stop_timer(bndry_flux_timer)

        end subroutine apply_bndry_fluxes

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_fill_halo
#ifndef ENABLE_DRY_MODE
            integer, parameter :: nfluxes = 2
#else
            integer, parameter :: nfluxes = 1
#endif
            call field_mpi_alloc(nfluxes, ndim=2)

            call field_interior_to_buffer(binc, 1)
#ifndef ENABLE_DRY_MODE
            call field_interior_to_buffer(hinc, 2)
#endif

            call interior_to_halo_communication

            call field_buffer_to_halo(binc, 1, .false.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_halo(hinc, 2, .false.)
#endif
            call field_mpi_dealloc

        end subroutine bndry_fluxes_fill_halo

end module bndry_fluxes
