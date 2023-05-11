!==============================================================================
! This module is used to apply boundary surface fluxes of buoyancy and
! humidity to each parcel in the lowest grid cell.
!==============================================================================
module bndry_fluxes
    use constants, only : zero, two
    use parameters, only : dx, lower, dxi, nx, ny
    use parcel_interpl, only : bilinear
    use fields, only : get_horizontal_index
    use parcel_container, only : n_parcels, parcels
    use netcdf_reader
    use omp_lib
    use mpi_utils, only : mpi_stop
    use field_mpi, only : field_mpi_alloc                   &
                        , field_mpi_dealloc                 &
                        , field_buffer_to_halo              &
                        , field_interior_to_buffer          &
                        , interior_to_halo_communication
    use mpi_layout, only : box
    implicit none

    logical, protected :: l_enable_flux
    logical            :: l_bndry_flux_allocated = .false.

    private

    ! Spatial form of the buoyancy and humidity fluxes through lower surface:
    double precision, dimension(:, :), allocatable :: bflux
#ifndef ENABLE_DRY_MODE
    double precision, dimension(:, :), allocatable :: hflux
#endif

    double precision :: fluxfac

    ! Denotes height below which surface fluxes are applied:
    double precision :: zdepth

    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! bilinear interpolation indices
    integer :: is(ngp), js(ngp)

    ! bilinear interpolation weights
    double precision :: weights(ngp)

    public :: l_enable_flux             &
            , apply_bndry_fluxes        &
            , read_bndry_fluxes         &
            , bndry_fluxes_allocate     &
            , bndry_fluxes_deallocate

    contains

        subroutine bndry_fluxes_allocate
            if (l_bndry_flux_allocated) then
                return
            endif

            l_bndry_flux_allocated = .true.

            allocate(thetaflux(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
#ifndef ENABLE_DRY_MODE
            allocate(qvflux(box%hlo(2):box%hhi(2), box%hlo(1):box%hhi(1)))
#endif
        end subroutine bndry_fluxes_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_deallocate
            if (.not. l_bndry_flux_allocated) then
                return
            endif

            l_bndry_flux_allocated = .false.

            deallocate(thetaflux)
#ifndef ENABLE_DRY_MODE
            deallocate(qvflux)
#endif
        end subroutine bndry_fluxes_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_default

            call bndry_fluxes_allocate

            bflux = zero
#ifndef ENABLE_DRY_MODE
            hflux = zero
#endif
        end subroutine bndry_fluxes_default

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine read_bndry_fluxes(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, start(3), cnt(3)
            integer                  :: lo(3), hi(3)

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

            if (has_dataset(ncid, 'thetaflux')) then
                call read_netcdf_dataset(ncid,                  &
                                         'thetaflux',               &
                                         thetaflux(lo(2):hi(2),     &
                                               lo(1):hi(1)),    &
                                         start,                 &
                                         cnt)
            else
                call mpi_stop("No theta flux field 'thetaflux' found in file.")
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'qvflux')) then
                call read_netcdf_dataset(ncid,                  &
                                         'qvflux',               &
                                         qvflux(lo(2):hi(2),     &
                                               lo(1):hi(1)),    &
                                         start,                 &
                                         cnt)
            else
                call mpi_stop("No qv flux field 'qvflux' found in file.")
            endif
#endif

            call close_netcdf_file(ncid)

            fluxfac = two * dxi(3) ** 2

            call bndry_fluxes_fill_halo

            !$omp parallel
            !$omp workshare
            thetaflux = fluxfac * thetaflux
#ifndef ENABLE_DRY_MODE
            qvflux = fluxfac * qvflux
#endif
            !$omp end workshare
            !$omp end parallel

            zdepth = lower(3) + dx(3)

        end subroutine read_bndry_fluxes

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Add fluxes of buoyancy and humidity to parcels in lowest grid layer:
        subroutine apply_bndry_fluxes(dt)
            double precision, intent(in)  :: dt
            double precision              :: xy(2), z, fac
            integer                       :: n, l

            if (.not. l_enable_flux) then
                return
            endif

            !$omp parallel default(shared)
            !$omp do private(n, l, is, js, weights, xy, z, fac)
            do n = 1, n_parcels
                z = parcels%position(3, n)
                if (z < zdepth) then

                    xy = parcels%position(1:2, n)

                    call bilinear(xy, is, js, weights)

                    fac = 2.0*(zdepth - z)/(dx(3)*dx(3)) * dt

                    do l = 1, ngp
                        ! The multiplication by dt is necessary to provide the amount of b or h
                        ! entering through the bottom surface over a time interval of dt.
                        parcels%theta(n) = parcels%theta(n) + fac * weights(l) * thetaflux(js(l), is(l))
#ifndef ENABLE_DRY_MODE
                        parcels%qv(n) = parcels%qv(n) + fac * weights(l) * qvflux(js(l), is(l))
#endif
                    enddo
                endif
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_bndry_fluxes

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine bndry_fluxes_fill_halo
#ifndef ENABLE_DRY_MODE
            integer, parameter :: nfluxes = 2
#else
            integer, parameter :: nfluxes = 1
#endif
            call field_mpi_alloc(nfluxes, ndim=2)

            call field_interior_to_buffer(thetaflux, 1)
#ifndef ENABLE_DRY_MODE
            call field_interior_to_buffer(qvflux, 2)
#endif

            call interior_to_halo_communication

            call field_buffer_to_halo(thetaflux, 1, .false.)
#ifndef ENABLE_DRY_MODE
            call field_buffer_to_halo(qvflux, 2, .false.)
#endif
            call field_mpi_dealloc

        end subroutine bndry_fluxes_fill_halo

end module bndry_fluxes
