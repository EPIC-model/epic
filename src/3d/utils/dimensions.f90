module dimensions
    implicit none

    ! number of dimensions
    integer, parameter :: n_dim = 3

    ! dimensional indices
    integer, parameter :: I_X    = 1 & ! index for zonal direction
                        , I_Y    = 2 & ! index for meridional direction
                        , I_Z    = 3   ! index for vertical direction

end module dimensions
