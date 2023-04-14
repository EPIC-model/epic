!==============================================================================
! This module is used to apply boundary surface fluxes of buoyancy and
! humidity to each parcel in the lowest grid cell.
!==============================================================================
module bndry_fluxes
    use constants, only : zero
    use parameters, only : nhcelli, dx, lower, dxi
    use parcel_interpl, only : bilinear
    use fields, only : get_index
    use parcel_container, only : n_parcels, parcels
    use omp_lib
    implicit none

    logical :: l_enable_flux = .false.

    private

    ! Spatial form of the buoyancy and humidity fluxes through lower surface:
    double precision, dimension(:, :), allocatable :: bflux, hflux

    ! Spatial average of the above quantities:
    double precision  :: avgbflux, avghflux

    double precision :: fluxfac

    ! Denotes height below which surface fluxes are applied:
    double precision :: zdepth

    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! bilinear interpolation indices
    integer :: is(ngp), js(ngp)

    ! bilinear interpolation weights
    double precision :: weights(ngp)

    contains

        subroutine read_bndry_fluxes
            if (.not. l_enable_flux) then
                return
            endif

            ! Compute spatial average:
            avgbflux = nhcelli * sum(bflux)

            ! Compute spatial average:
            avghflux = nhcelli * sum(hflux)

            fluxfac = two * dxi(3) ** 2

            !$omp parallel
            !$omp workshare
            bflux = fluxfac * bflux
            hflux = fluxfac * hflux
            !$omp end workshare
            !$omp end parallel

            zdepth = lower(3) + dx(3)

        end subroutine read_bndry_fluxes

        ! Add fluxes of buoyancy and humidity to parcels in lowest grid layer:
        subroutine apply_bndry_fluxes(dt)
            double precision, intent(in)  :: dt
            double precision              :: xy(2), z
            double precision, allocatable :: bf(:)
            double precision              :: dbm, bfsum
#ifndef ENABLE_DRY_MODE
            double precision, allocatable :: hf(:)
            double precision              :: dhm, hfsum
#endif
            integer                       :: n, i, j, l

            if (.not. l_enable_flux) then
                return
            endif

            bfsum = zero
#ifndef ENABLE_DRY_MODE
            hfsum = zero
#endif

#ifndef ENABLE_DRY_MODE
            !$omp parallel default(shared)
            !$omp do private(n, i, j, l, is, js, weights, xy, z) &
            !$omp& reduction(+: bfsum, hfsum)
#else
            !$omp parallel default(shared)
            !$omp do private(n, i, j, l, is, js, weights, xy, z) &
            !$omp& reduction(+: bfsum)
#endif
            do n = 1, n_parcels
                z = parcels%position(3, n)
                if (z < zdepth) then

                    xy = parcels%position(1:2, n)

                    call get_horizontal_index(xy, i, j)

                    call bilinear(xy, is, js, weights)

                    do l = 1, ngp
                        bf(n) = bf(n) + weights(l) * bflux(js(l), is(l))
#ifndef ENABLE_DRY_MODE
                        hf(n) = hf(n) + weights(l) * hflux(js(l), is(l))
#endif
                    enddo

                    ! buoynacy flux increment:
                    bf(n) = (zdepth - z) * bf(n)
                    bfsum = bfsum + bf(n) * parcels%volume(n)

#ifndef ENABLE_DRY_MODE
                    ! humidity flux increment:
                    hf(n) = (zdepth - z) * hf(n)
                    hfsum = hfsum + hf(n) * parcels%volume(n)
#endif
                endif
            enddo
            !$omp end do
            !$omp end parallel

            ! Complete calculation of average fluxes on parcels:
            bfsum = bfsum * dx(3) * nhcelli
#ifndef ENABLE_DRY_MODE
            hfsum = hfsum * dx(3) * nhcelli
#endif

            ! Correct parcels fluxes to be consistent with prescribed spatial fluxes:
            ! The multiplication by dt is necessary to provide the amount of b or h
            ! entering through the bottom surface over a time interval of dt.
            dbm = dt * avgbflux / bfsum
#ifndef ENABLE_DRY_MODE
            dhm = dt * avghflux / hfsum
#endif

            !$omp parallel default(shared)
            !$omp do private(n)
            do n = 1, n_parcels
                if (parcels%position(3, n) < zdepth) then
                    ! Scale b & h values to give correct net flux:
                    parcels%buoyancy(n) = parcels%buoyancy(n) + dbm * bf(n)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(n) = parcels%humidity(n) + dhm * hf(n)
#endif
                endif
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_bndry_fluxes

end module bndry_fluxes
