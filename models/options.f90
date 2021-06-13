module options
    implicit none

    !
    ! TaylorGreen options
    !


    !
    ! Straka case options
    !


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
