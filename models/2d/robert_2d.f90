! =============================================================================
!                       Robert ("Rising bubble") test case
!
!   This module sets up a bubble with uniform or Gaussian buoyancy profile.
!   The setup is from the paper
!
!   Robert, A. (1993). Bubble Convection Experiments with a Semi-implicit
!   Formulation of the Euler Equations, Journal of Atmospheric Sciences,
!   50(13), 1865-1873. Retrieved Jun 2, 2021, from
!   https://doi.org/10.1175/1520-0469(1993)050<1865:BCEWAS>2.0.CO;2
!
! =============================================================================

module robert_2d
    use constants
    use netcdf_writer
    use physical_constants, only : write_physical_constants, gravity
    use physical_parameters, only : write_physical_parameters, theta_0
    implicit none

    private

    double precision, allocatable :: buoyg(:, :)

    integer :: buo_id

    type bubble_type
        character(len=8) :: distr           ! distribution ('gaussian' or 'uniform')
        double precision :: center(2)       ![m] bubble center (x, z)
        double precision :: dtheta_max      ![K] max. pot. temp. perturbation
        double precision :: radius          ![m] bubble outer radius ('uniform')
                                            ! or plateau radius ('gaussian')
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

        subroutine robert_init(ncid, dimids, nx, nz, origin, dx)
            integer,          intent(inout) :: ncid
            integer,          intent(in)    :: dimids(:)
            integer,          intent(in)    :: nx, nz
            double precision, intent(in)    :: origin(2)
            double precision, intent(in)    :: dx(2)
            integer                         :: k
            type(bubble_type)               :: bubble

            ! set physical parameters
            theta_0 = 303.15d0         ![K] reference potential temperature

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

            allocate(buoyg(0:nz, 0:nx-1))

            buoyg = zero

            do k = 1, robert_flow%n_bubbles
                bubble = robert_flow%bubbles(k)

                select case (trim(bubble%distr))
                    case ('uniform')
                        call robert_uniform_init(nx, nz, origin, dx, bubble)
                    case ('gaussian')
                        call robert_gaussian_init(nx, nz, origin, dx, bubble)
                    case default
                        print *, "Invalid distribution type: '", trim(bubble%distr), "'"
                        stop
                end select
            enddo

            call write_netcdf_dataset(ncid, buo_id, buoyg)

            call write_physical_constants(ncid)
            call write_physical_parameters(ncid)

            deallocate(buoyg)

        end subroutine robert_init

        subroutine robert_uniform_init(nx, nz, origin, dx, bubble)
            integer,           intent(in) :: nx, nz
            double precision,  intent(in) :: origin(2), dx(2)
            type(bubble_type), intent(in) :: bubble
            double precision              :: xc, zc, r2, r2_out, dtheta, pos(2)
            integer                       :: i, j

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r2_out = bubble%radius ** 2

            do j = 0, nz
                do i = 0, nx-1
                    pos = origin + dx * dble((/i, j/))

                    r2 = (pos(1) - xc) ** 2 &
                       + (pos(2) - zc) ** 2

                    ! potential temperature perturbation
                    dtheta = zero

                    if (r2 <= r2_out) then
                        dtheta = bubble%dtheta_max
                    endif

                    ! MPIC paper:
                    ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                    ! (dtheta = theta - theta_0)
                    buoyg(j, i) = buoyg(j, i) &
                                + gravity * dtheta / theta_0
                enddo
            enddo
        end subroutine robert_uniform_init


        subroutine robert_gaussian_init(nx, nz, origin, dx, bubble)
            integer,           intent(in) :: nx, nz
            double precision,  intent(in) :: origin(2), dx(2)
            type(bubble_type), intent(in) :: bubble
            integer                       :: i, j
            double precision              :: xc, zc, r, r_in, dtheta
            double precision              :: dtheta_max, pos(2), s, fs2

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r_in = bubble%radius
            s = bubble%width
            fs2 = (one / s) ** 2

            ! max. potential temperature perturbation
            dtheta_max = bubble%dtheta_max

            do j = 0, nz
                do i = 0, nx-1
                    pos = origin + dx * dble((/i, j/))

                    r = (pos(1) - xc) ** 2 &
                      + (pos(2) - zc) ** 2

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
                    buoyg(j, i) = buoyg(j, i) &
                                + gravity * dtheta / theta_0
                enddo
            enddo
        end subroutine robert_gaussian_init
end module
