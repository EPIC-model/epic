module constants
    implicit none

    double precision, parameter :: zero = 0.d0

    double precision, parameter :: one = 1.d0

    double precision, parameter :: two = 2.d0

    double precision, parameter :: pi = acos(-one)

    ! maximum number of allowed parcels
    integer, parameter :: max_num_parcels = 1e6

    ! domain size
    double precision, parameter :: extent(2) = (/ pi, pi /)

    ! domain half widths and edge values:
    double precision, parameter :: hl(2) = extent / two

    double precision, parameter :: hli(2) = one / hl

    ! domain origin
    double precision, parameter :: lower(2) = - hl

    ! domain upper boundary
    double precision, parameter :: upper(2) = hl

end module
