module macros
    implicit none

#ifdef ENABLE_3D
    integer, parameter   :: EPIC_DIM = 3
    character, parameter :: EPIC_DIM_CHAR = '3'
#else
    integer, parameter   :: EPIC_DIM = 2
    character, parameter :: EPIC_DIM_CHAR = '2'
#endif

    interface EPIC_VECTOR
        module procedure :: EPIC_DVECTOR
        module procedure :: EPIC_IVECTOR
    end interface EPIC_VECTOR

    contains

#ifdef ENABLE_3D
    !
    !   3D version
    !
    pure function EPIC_DVECTOR(i, j, k) result(vec)
        double precision, intent(in) :: i, j, k
        double precision             :: vec(EPIC_DIM)

        vec = (/i, j, k/)
    end function EPIC_DVECTOR

    pure function EPIC_IVECTOR(i, j, k) result(vec)
        integer, intent(in) :: i, j, k
        integer             :: vec(EPIC_DIM)

        vec = (/i, j, k/)
    end function EPIC_IVECTOR

#else
    !
    !   2D version
    !
    pure function EPIC_DVECTOR(i, j, k) result(vec)
        double precision,           intent(in) :: i, j
        double precision, optional, intent(in) :: k
        double precision                       :: vec(EPIC_DIM)

        vec = (/i, j/)
    end function EPIC_DVECTOR

    pure function EPIC_IVECTOR(i, j, k) result(vec)
        integer,           intent(in) :: i, j
        integer, optional, intent(in) :: k
        integer                       :: vec(EPIC_DIM)

        vec = (/i, j/)
    end function EPIC_IVECTOR
#endif
end module macros
