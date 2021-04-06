module constants
    implicit none

    double precision, parameter :: zero = 0.d0

    double precision, parameter :: one = 1.d0

    double precision, parameter :: two = 2.d0

    double precision, parameter :: pi = acos(-one)

    ! maximum number of allowed parcels
    integer, parameter :: max_num_parcels = 1e6

end module
