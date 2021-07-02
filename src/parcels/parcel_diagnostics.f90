! =============================================================================
!                               Parcel diagnostics
! =============================================================================
module parcel_diagnostics
    use constants, only : zero, one, f12
    use merge_sort
    use parameters, only : extent, lower
    use parcel_container, only : parcels, n_parcels
    use h5_utils
    use h5_writer
    use omp_lib
    implicit none

    private

    ! peref : potential energy reference
    ! pe    : potential energy
    ! ke    : kinetic energy
    double precision :: peref, pe, ke

    public :: init_parcel_diagnostics, &
              write_h5_parcel_diagnostics

    contains

        ! compute the reference potential energy
        subroutine init_parcel_diagnostics
            integer          :: ii(n_parcels), n
            double precision :: b(n_parcels)
            double precision :: gam, zmean

            b = parcels%buoyancy(1:n_parcels)

            ! sort buoyancy in ascending order
            call msort(b, ii)

            gam = one / extent(1)
            zmean = lower(2) + f12 * gam * parcels%volume(ii(1))

            peref = - b(1) * parcels%volume(ii(1)) * zmean
            do n = 2, n_parcels
                zmean = zmean + gam * parcels%volume(ii(n))

                peref = peref &
                      - b(n) * parcels%volume(ii(n)) * zmean
            enddo
        end subroutine init_parcel_diagnostics


        subroutine calculate_diagnostics
            integer          :: n
            double precision :: b, z, v, vel(2), vol, zmin

            ! reset
            ke = zero
            pe = zero

            zmin = lower(2)

            !$omp parallel
            !$omp do private(n, vel, vol, b, z) reduction(+: ke, pe)
            do n = 1, n_parcels

                vel = parcels%velocity(n, :)
                vol = parcels%volume(n)
                b   = parcels%buoyancy(n)
                z   = parcels%position(n, 2) - zmin

                ! kinetic energy
                ke = ke + (vel(1) ** 2 + vel(2) ** 2) * vol

                ! potential energy
                pe = pe - b * z * vol
            enddo
            !$omp end do
            !$omp end parallel

            ke = f12 * ke
            pe = pe - peref
        end subroutine calculate_diagnostics


        subroutine write_h5_parcel_diagnostics(h5file_id, iter)
            integer(hid_t), intent(in)    :: h5file_id
            integer,        intent(in)    :: iter ! iteration
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            logical                       :: created

            call calculate_diagnostics
            name = trim(get_step_group_name(iter))

            call create_h5_group(h5file_id, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id, name, group)
            endif

            !
            ! write diagnostics
            !
            call write_h5_double_scalar_attrib(group, "potential energy", pe)
            call write_h5_double_scalar_attrib(group, "kinetic energy", ke)
            call write_h5_double_scalar_attrib(group, "total energy", ke + pe)

            ! close all
            call close_h5_group(group)
        end subroutine write_h5_parcel_diagnostics
end module parcel_diagnostics
