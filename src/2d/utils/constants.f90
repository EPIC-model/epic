! =============================================================================
! This module contains all non-modifiable parameters and all quantities which
! never change throughout a simulation.
! =============================================================================
module constants
    implicit none

    double precision, parameter :: zero     = 0.0d0
    double precision, parameter :: one      = 1.0d0
    double precision, parameter :: two      = 2.0d0
    double precision, parameter :: three    = 3.0d0
    double precision, parameter :: four     = 4.0d0
    double precision, parameter :: five     = 5.0d0
    double precision, parameter :: six      = 6.0d0
    double precision, parameter :: ten      = 10.0d0
    double precision, parameter :: hundred  = ten ** 2
    double precision, parameter :: thousand = ten ** 3

    double precision, parameter :: pi    = dacos(-one)
    double precision, parameter :: twopi = two * pi
    double precision, parameter :: fpi   = one / pi

    double precision, parameter :: f32   = three / two
    double precision, parameter :: f34   = three / four
    double precision, parameter :: f12   = one / two
    double precision, parameter :: f13   = one / three
    double precision, parameter :: f14   = one / four
    double precision, parameter :: f23   = two / three
    double precision, parameter :: f112  = one / 12.d0
    double precision, parameter :: f16   = one / six
    double precision, parameter :: f56   = five / six
    double precision, parameter :: f76   = 7.d0 / six
    double precision, parameter :: f124  = one / 24.d0
    double precision, parameter :: f1112 = 11.d0 / 12.d0

    ! maximum number of allowed parcels
    integer, parameter :: max_num_parcels = 2.2e6

end module
