module unit_test
    implicit none

    ! abolute tolerance
    double precision, parameter :: atol_m = 1.0e-15

    private :: atol_m

    contains


        ! test : test name
        ! aerr : absolute error
        ! atol : absolute tolerance (optional, default: 1e-15)
        subroutine print_result_dp(test, aerr, atol)
            character(*),               intent(in) :: test
            double precision,           intent(in) :: aerr
            double precision, optional, intent(in) :: atol
            double precision                       :: tol

            tol = atol_m

            if (present(atol)) then
                tol = atol
            endif

            call print_result_logical(test, aerr > tol)

        end subroutine print_result_dp

        subroutine print_result_logical(test, failed)
            character(*),    intent(in) :: test
            logical,         intent(in) :: failed
            character(len=6)            :: res

            res = 'PASSED'
            if (failed) then
                res = 'FAILED'
            endif
            ! 5 Mai 2021
            ! https://stackoverflow.com/questions/47761900/format-add-trailing-spaces-to-character-output-to-left-justify
            print '(a60, a7)', ' '//test//':'//repeat(' ', 57), res

        end subroutine print_result_logical
end module unit_test
