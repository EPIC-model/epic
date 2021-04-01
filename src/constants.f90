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

    ! domain origin
    double precision, parameter :: lower(2) = - extent / two

    ! domain upper boundary
    double precision, parameter :: upper(2) = extent / two
end module
