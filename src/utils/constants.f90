! =============================================================================
! This module contains all non-modifiable parameters and all quantities which
! never change throughout a simulation.
! =============================================================================
module constants
    implicit none

    double precision, parameter :: zero = 0.0d0

    double precision, parameter :: one = 1.0d0

    double precision, parameter :: two = 2.0d0

    double precision, parameter :: three = 3.0d0

    double precision, parameter :: four = 4.0d0

    double precision, parameter :: pi = acos(-one)

    double precision, parameter :: twopi = two * pi

    double precision, parameter :: one_over_pi = one / pi

    ! maximum number of allowed parcels
    integer, parameter :: max_num_parcels = 1e6

end module
