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

module robert
    use physics
    use constants
    use writer
    implicit none

    private

    double precision, allocatable :: buoyg(:, :)

    type bubble_type
        character(len=8) :: distr           ! distribution ('gaussian' or 'uniform')
        double precision :: center(2)       ![m] bubble center (x, z)
        double precision :: theta_max       ![Kelvin] max. pot. temp. perturbation
        double precision :: outer_radius    ![m] bubble outer radius
        double precision :: inner_radius    ![m] bubble inner radius (if 'gaussian')
        double precision :: width           ![m] standard deviation of Gaussian
    end type bubble_type

    type flow_type
        double precision  :: theta_ref   = 303.15d0   ![Kelvin] reference potential temperature
        integer           :: n_bubbles   = 1
        type(bubble_type) :: bubbles(10)
    end type flow_type

    type(flow_type) :: robert_flow

    public :: robert_init, &
              robert_flow

    contains

        subroutine robert_init(filename, nx, nz, origin, dx)
            character(*),     intent(in) :: filename
            integer,          intent(in) :: nx, nz
            double precision, intent(in) :: origin(2)
            double precision, intent(in) :: dx(2)
            integer           :: k
            type(bubble_type) :: bubble

            allocate(buoyg(0:nz, 0:nx-1))

            if (robert_flow%n_bubbles > size(robert_flow%bubbles)) then
                print *, 'Number of bubbles beyond upper limit.'
                stop
            endif


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

            call open_h5_file(filename)
            call write_h5_dataset_2d('/', 'buoyancy', buoyg)
            call close_h5_file

            deallocate(buoyg)

        end subroutine robert_init

        subroutine robert_uniform_init(nx, nz, origin, dx, bubble)
            integer,           intent(in) :: nx, nz
            double precision,  intent(in) :: origin(2), dx(2)
            type(bubble_type), intent(in) :: bubble
            double precision              :: xc, zc, r2, r2_out, dtheta, theta_ref, pos(2)
            integer                       :: i, j

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r2_out = bubble%outer_radius ** 2

            ! reference potential temperature
            theta_ref = robert_flow%theta_ref

            do j = 0, nz
                do i = 0, nx-1
                    pos = origin + dx * dble((/i, j/))

                    r2 = (pos(1) - xc) ** 2 &
                       + (pos(2) - zc) ** 2

                    ! potential temperature perturbation
                    dtheta = zero

                    if (r2 <= r2_out) then
                        dtheta = bubble%theta_max
                    endif

                    ! MPIC paper:
                    ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                    ! (dtheta = theta - theta_ref)
                    buoyg(j, i) = gravity * dtheta / theta_ref
                enddo
            enddo
        end subroutine robert_uniform_init


        subroutine robert_gaussian_init(nx, nz, origin, dx, bubble)
            integer,           intent(in) :: nx, nz
            double precision,  intent(in) :: origin(2), dx(2)
            type(bubble_type), intent(in) :: bubble
            integer                       :: i, j
            double precision              :: xc, zc, r, r_out, r_in, dtheta
            double precision              :: theta_ref, theta_max, pos(2), s, fs2

            ! in metres
            xc = bubble%center(1)
            zc = bubble%center(2)
            r_out = bubble%outer_radius
            r_in = bubble%inner_radius
            s = bubble%width
            fs2 = (one / s) ** 2

            ! reference potential temperature
            theta_ref = robert_flow%theta_ref

            ! max. potential temperature perturbation
            theta_max = bubble%theta_max

            do j = 0, nz
                do i = 0, nx-1
                    pos = origin + dx * dble((/i, j/))

                    r = (pos(1) - xc) ** 2 &
                      + (pos(2) - zc) ** 2

                    r = dsqrt(r)

                    ! potential temperature perturbation
                    dtheta = zero

                    if (r <= r_out) then
                        if (r <= r_in) then
                            dtheta = theta_max
                        else
                            dtheta = theta_max * dexp(-(r - r_in) ** 2 * fs2)
                        endif
                    endif

                    ! MPIC paper:
                    ! liquid-water buoyancy is defined by b = g * (theta − theta_ref) / theta_ref
                    ! (dtheta = theta - theta_ref)
                    buoyg(j, i) = gravity * dtheta / theta_ref
                enddo
            enddo
        end subroutine robert_gaussian_init
end module
