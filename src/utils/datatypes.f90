module datatypes
    ! See also https://fortran-lang.org/en/learn/intrinsics/type/#selected-int-kind
    ! for reference
    implicit none

    ! 8-byte integer
    ! Use "integer(kind=int64)" to create an integer of size 8 bytes
    ! The maximum value of a signed integer is 2^63-1 = 9,223,372,036,854,775,807,
    ! i.e. we must be able to store values from -10^18 to 10^18:
    integer, parameter :: int64 = selected_int_kind(18)

    type :: intlog_pair_t
        sequence
        integer :: ival
        logical :: lval
    end type intlog_pair_t

end module datatypes
