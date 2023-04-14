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
    implicit none

    logical, protected :: l_enable_flux

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

    public :: l_enable_flux, apply_bndry_fluxes, read_bndry_fluxes

    contains

        subroutine read_bndry_fluxes(fname)
            character(*), intent(in) :: fname
            integer                  :: ncid, start(3), cnt(3)

            l_enable_flux = (fname /= '')

            if (.not. l_enable_flux) then
                return
            endif

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            cnt  =  (/ nx, ny, 1 /)
            start = (/ 1,  1,  1 /)

            if (has_dataset(ncid, 'bflux')) then
                allocate(bflux(ny, nx))
                call read_netcdf_dataset(ncid, 'bflux', bflux, start, cnt)
            else
                print *, "No buoyancy flux field 'bflux' found in file."
                stop
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'hflux')) then
                allocate(hflux(ny, nx))
                call read_netcdf_dataset(ncid, 'hflux', hflux, start, cnt)
            else
                print *, "No humidity flux field 'hflux' found in file."
                stop
            endif
#endif

            call close_netcdf_file(ncid)

            fluxfac = two * dxi(3) ** 2

            !$omp parallel
            !$omp workshare
            bflux = fluxfac * bflux
#ifndef ENABLE_DRY_MODE
            hflux = fluxfac * hflux
#endif
            !$omp end workshare
            !$omp end parallel

            zdepth = lower(3) + dx(3)

        end subroutine read_bndry_fluxes

        ! Add fluxes of buoyancy and humidity to parcels in lowest grid layer:
        subroutine apply_bndry_fluxes(dt)
            double precision, intent(in)  :: dt
            double precision              :: xy(2), z, fac
            integer                       :: n, i, j, l

            if (.not. l_enable_flux) then
                return
            endif

            !$omp parallel default(shared)
            !$omp do private(n, i, j, l, is, js, weights, xy, z)
            do n = 1, n_parcels
                z = parcels%position(3, n)
                if (z < zdepth) then

                    xy = parcels%position(1:2, n)

                    call get_horizontal_index(xy, i, j)

                    call bilinear(xy, is, js, weights)

                    fac = (zdepth - z) * dt

                    do l = 1, ngp
                        ! The multiplication by dt is necessary to provide the amount of b or h
                        ! entering through the bottom surface over a time interval of dt.
                        parcels%buoyancy(n) = parcels%buoyancy(n) + fac * weights(l) * bflux(js(l), is(l))
#ifndef ENABLE_DRY_MODE
                        parcels%humidity(n) = parcels%humidity(n) + fac * weights(l) * hflux(js(l), is(l))
#endif
                    enddo
                endif
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_bndry_fluxes

end module bndry_fluxes
