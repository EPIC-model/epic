module options
    implicit none

    !
    ! TaylorGreen options
    !
    type taylor_green_type
        double precision :: amp(2) = 1.0d0    ! amplitudes
        double precision :: freq(2) = 1.0d0   ! frequencies
        double precision :: phase(2) = 1.0d0  ! phase shift
    end type taylor_green_type

    type(taylor_green_type) :: taylor_green_opt

    !
    ! Straka case options
    !
    type straka_type
        double precision :: theta_ref = 300.0d0               ![Kelvin] reference potential temperature
        double precision :: theta_max = 15.0d0                ![Kelvin] max. pot. temp. perturbation
        double precision :: center(2) = (/0.0d0, 3000.0d0/)   ![m] sphere center (x, z)
        double precision :: radii(2)  = (/4000.0d0, 2000.d0/) ![m] ellipse radii (x, z)
    end type straka_type

    type(straka_type) :: straka_opt

    !
    ! Robert case options
    !
    type bubble_type
        character(len=8) :: distr           ! distribution ('gaussian' or 'uniform')
        double precision :: center(2)       ![m] bubble center (x, z)
        double precision :: theta_max       ![Kelvin] max. pot. temp. perturbation
        double precision :: outer_radius    ![m] bubble outer radius
        double precision :: inner_radius    ![m] bubble inner radius (if 'gaussian')
        double precision :: width           ![m] standard deviation of Gaussian
    end type bubble_type

    type robert_type
        double precision  :: theta_ref   = 303.15d0   ![Kelvin] reference potential temperature
        integer           :: n_bubbles   = 1
        type(bubble_type) :: bubbles(10)
    end type robert_type

    type(robert_type) :: robert_opt

end module options
