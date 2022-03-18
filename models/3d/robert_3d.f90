! =============================================================================
!                       Robert ("Rising bubble") test case
!
!   This module sets up to ten bubbles with a flat-top ("plateau") beyond
!   which the buoyancy follows a Gaussian profile.
! =============================================================================

module robert_3d
    use constants
    use netcdf_writer
    use physics, only : write_physical_quantities, &
                        set_physical_quantity, gravity, theta_0
    implicit none

    private

    double precision, allocatable :: buoyg(:, :, :)

    integer :: buo_id

    type bubble_type
        double precision :: center(3)       ![m] bubble center (x, z)
        double precision :: dtheta_max      ![K] max. pot. temp. perturbation
        double precision :: radius          ![m] plateau radius
        double precision :: width           ![m] standard deviation of Gaussian
    end type bubble_type

    type flow_type
        integer           :: n_bubbles   = 1
        type(bubble_type) :: bubbles(10)
    end type flow_type

    type(flow_type) :: robert_flow

    public :: robert_init, &
              robert_flow

    contains

        subroutine robert_init(ncid, dimids, nx, ny, nz, origin, dx)
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(in)    :: origin(3)
            double precision, intent(in)    :: dx(3)
            integer                         :: k
            type(bubble_type)               :: bubble

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='buoyancy',                     &
                                       long_name='buoyancy',                &
                                       std_name='',                         &
                                       unit='m/s^2',                        &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=buo_id)

            call close_definition(ncid)

            if (robert_flow%n_bubbles > size(robert_flow%bubbles)) then
                print *, 'Number of bubbles beyond upper limit.'
                stop
            endif

            allocate(buoyg(0:nz, 0:ny-1, 0:nx-1))

            buoyg = zero

            do k = 1, robert_flow%n_bubbles
                bubble = robert_flow%bubbles(k)

                call robert_gaussian_init(nx, ny, nz, origin, dx, bubble)
            enddo

            call write_netcdf_dataset(ncid, buo_id, buoyg)

            call write_physical_quantities(ncid)

            deallocate(buoyg)

        end subroutine robert_init

        subroutine robert_gaussian_init(nx, ny, nz, origin, dx, bubble)
            integer,           intent(in) :: nx, ny, nz
            double precision,  intent(in) :: origin(3), dx(3)
            type(bubble_type), intent(in) :: bubble
            integer                       :: i, j, k
            double precision              :: xc, yc, zc, r, r_in, dtheta
            double precision              :: dtheta_max, pos(3), s, fs2

            ! in metres
            xc = bubble%center(1)
            yc = bubble%center(2)
            zc = bubble%center(3)
            r_in = bubble%radius
            s = bubble%width
            fs2 = (one / s) ** 2

            ! max. potential temperature perturbation
            dtheta_max = bubble%dtheta_max

            do i = 0, nx-1
                do j = 0, ny-1
                    do k = 0, nz
                        pos = origin + dx * dble((/i, j, k/))

                        r = (pos(1) - xc) ** 2 &
                          + (pos(2) - yc) ** 2 &
                          + (pos(3) - zc) ** 2

                        r = dsqrt(r)

                        ! potential temperature perturbation
                        dtheta = zero

                        if (r <= r_in) then
                            dtheta = dtheta_max
                        else
                            dtheta = dtheta_max * dexp(-(r - r_in) ** 2 * fs2)
                        endif

                        ! MPIC paper:
                        ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                        ! (dtheta = theta - theta_0)
                        buoyg(k, j, i) = buoyg(k, j, i) &
                                       + gravity * dtheta / theta_0
                    enddo
                enddo
            enddo
        end subroutine robert_gaussian_init
end module robert_3d
