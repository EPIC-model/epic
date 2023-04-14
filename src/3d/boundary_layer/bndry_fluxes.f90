!==============================================================================
! This module is used to apply boundary surface fluxes of buoyancy and
! humidity to each parcel in the lowest grid cell.
!==============================================================================
module bndry_fluxes
    use constants, only : zero
    use parameters, only : nhcelli, dx
    use parcels_interpl, only : bilinear
    use fields, only : get_index
    implicit none

    logical :: l_enable_flux = .false.

    private

    ! Spatial form of the buoyancy and humidity fluxes through lower surface:
    double precision, dimension(:, :), allocatable :: bflux, hflux

    ! Spatial average of the above quantities:
    double precision  :: avgbflux, avghflux

    double precision :: fluxfac! = two*dzi**2

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

            !$omp parallel
            !$omp workshare
            bflux = fluxfac * bflux
            hflux = fluxfac * hflux
            !$omp end workshare
            !$omp end parallel

            zdepth = lower(3) + dx(3)


        end subroutine read_bndry_fluxes

        subroutine apply_bndry_fluxes(dt)
            double precision, intent(in)  :: dt
            double precision              :: xy(2), z, dbm, dhm
            double precision, allocatable :: bf(:), hf(:)

            if (.not. l_enable_flux) then
                return
            endif

            ! Add fluxes of buoyancy and humidity to parcels in lowest grid layer:
            bfsum = zero
            hfsum = zero

            !$omp parallel default(shared)
            !$omp do private(n, i, j, l, is, js, weights, xy, z) schedule(static, chunk) &
            !$omp& reduction(+: bfsum, hfsum)
            do n = 1, n_parcels
                z = parcels%position(3, n)
                if (z < zdepth) then

                    xy = parcels%position(1:2, n)

                    call get_index(xy, i, j)

                    call bilinear(xy, is, js weights)

                    do l = 1, ngp
                        bf(n) = bf(n) + weights(l) * bflux(js(l), is(l))
                        hf(n) = hf(n) + weights(l) * hflux(js(l), is(l))
                    enddo

                    ! buoynacy flux increment:
                    bf(n) = (zdepth - z) * bf(n)
                    bfsum = bfsum + bf(n) * parcels%volume(n)

                    ! humidity flux increment:
                    hf(n) = (zdepth - z) * hf(n)
                    hfsum = hfsum + hf(n) * parcels%volume(n)
                endif
            enddo
            !$omp end do
            !$omp end parallel

            ! Complete calculation of average fluxes on parcels:
            bfsum = bfsum * dx(3) * nhcelli
            hfsum = hfsum * dx(3) * nhcelli

            ! Correct parcels fluxes to be consistent with prescribed spatial fluxes:
            ! The multiplication by dt is necessary to provide the amount of b or h
            ! entering through the bottom surface over a time interval of dt.
            dbm = dt * avgbflux / bfsum
            dhm = dt * avghflux / hfsum

            !$omp parallel default(shared)
            !$omp do private(n) schedule(static, chunk)
            do n = 1, n_parcels
                if (parcels%position(3, n) < zdepth) then
                    ! Scale b & h values to give correct net flux:
                    parcels%buoyancy(n) = parcels%buoyancy(n) + dbm * bf(n)
                    parcels%humidity(n) = parcels%humidity(n) + dhm * hf(n)
                endif
            enddo
            !$omp end do
            !$omp end parallel

        end subroutine apply_bndry_fluxes

end module bndry_fluxes
