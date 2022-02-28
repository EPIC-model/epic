module unit_test
    implicit none

    ! abolute tolerance
    double precision, parameter :: atol_m = dble(1.0e-15)

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

            call print_result_logical(test, aerr < tol)

        end subroutine print_result_dp

        subroutine print_result_logical(test, passed)
            character(*),    intent(in) :: test
            logical,         intent(in) :: passed
            character(len=17)           :: res

            ! 5 August 2021
            ! https://stackoverflow.com/questions/6402700/coloured-terminal-output-from-fortran
            ! https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
            res = achar(27) // '[32;1mPASSED' // achar(27) // '[0m'
            if (.not. passed) then
                res = achar(27) // '[31;1mFAILED' // achar(27) // '[0m'
            endif
            ! 5 Mai 2021
            ! https://stackoverflow.com/questions/47761900/format-add-trailing-spaces-to-character-output-to-left-justify
            print '(a60, a18)', ' '//test//' '//repeat('.', 57), res

        end subroutine print_result_logical

        ! Give option to run unit tests verbosely using command line
        subroutine parse_command_line
            integer           :: i
            character(len=32) :: arg
            character(len=32) :: testname
            logical           :: verbose

            testname = ''
            i = 0
            do
                call get_command_argument(i, arg)
                if(i==0) then
                    testname=trim(arg)
                elseif (len_trim(arg) == 0) then
                    exit
                elseif (arg == '--verbose') then
                    verbose = .true.
                endif
                i = i+1
            end do

            if (verbose) then
                print *, 'Running ', trim(testname),' verbosely'
            endif
        end subroutine parse_command_line

end module unit_test
