! =============================================================================
! This module contains data type definitions
! =============================================================================
module datatypes
    use mpi
    !use iso_fortran_env
    implicit none
	    
    integer, public, parameter :: STRING_LENGTH=150 !< Default length of strings
    integer, public, parameter :: LONG_STRING_LENGTH=STRING_LENGTH + 50!< Length of longer strings

    integer, public, parameter :: SINGLE_PRECISION = selected_real_kind(6,30)   ! Single precision (32 bit) kind
    integer, public, parameter :: DOUBLE_PRECISION = selected_real_kind(15,307) ! Double precision (64 bit) kind

    integer, public, parameter :: LONG_INTEGER_KIND = 8
    integer, public, parameter :: SHORT_INTEGER_KIND = 4

    integer, public, parameter :: PARCEL_INTEGER = MPI_INTEGER8 !LONG_INTEGER

    ! Default precision for prognostic and calculations data
    integer, public, parameter :: DEFAULT_PRECISION = DOUBLE_PRECISION
	
    ! MPI communication type for the prognostic and calculation data
    integer, public :: PRECISION_TYPE = MPI_DOUBLE_PRECISION
    integer, public :: MPI_PARCEL_INT = MPI_INTEGER8

end module
