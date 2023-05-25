module linalg
    use constants, only : zero, one
    implicit none

    contains

    ! sign function (may be put in another module)
    pure function signum(val) result(sgn)
        double precision, intent(in) :: val
        double precision             :: sgn

        sgn = sign(one, val)

        if (val == zero) then
            sgn = zero
        endif
    end function signum

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates the cross product of two vectors 'a' and 'b'.
    pure function cross(a, b) result(c)
        double precision, intent(in) :: a(3), b(3)
        double precision             :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates the trace of a matrix A.
    pure function trace(A) result(tr)
        double precision, intent(in) :: A(3, 3)
        double precision             :: tr

        tr = A(1, 1) + A(2, 2) + A(3, 3)

    end function trace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Calculates the determinant of a matrix A.
    pure function determinant(A) result(det)
        double precision, intent(in) :: A(3, 3)
        double precision             :: det

        det = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) ** 2)      &
            - A(1, 2) * (A(1, 2) * A(3, 3) - A(1, 3) * A(2, 3)) &
            + A(1, 3) * (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2))

    end function determinant

end module linalg
