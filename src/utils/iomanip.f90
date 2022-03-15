module iomanip

    interface print_key_value_pair
        module procedure :: print_key_value_pair_double
        module procedure :: print_key_value_pair_integer
        module procedure :: print_key_value_pair_logical
        module procedure :: print_key_value_pair_character
    end interface print_key_value_pair

    private :: print_key_value_pair_double,     &
               print_key_value_pair_integer,    &
               print_key_value_pair_logical,    &
               print_key_value_pair_character
    contains

    ! convert number to string of length 10 with
    ! leading zeros
    function zfill(num) result(name)
        integer, intent(in) :: num
        ! 12 March 2021
        ! https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
        character(len=10) :: name

        write(name, fmt='(I10.10)') num
    end function zfill


    subroutine print_key_value_pair_double(name, val)
        character(*),     intent(in) :: name
        double precision, intent(in) :: val
        character(40)                :: fix_length_name

        fix_length_name = name
        write (*, "(a, '=', f15.3)") fix_length_name, val
    end subroutine print_key_value_pair_double

    subroutine print_key_value_pair_integer(name, val)
        character(*), intent(in) :: name
        integer,      intent(in) :: val
        character(40)            :: fix_length_name

        fix_length_name = name
        write (*, "(a, '=', I15)") fix_length_name, val
    end subroutine print_key_value_pair_integer

    subroutine print_key_value_pair_logical(name, val)
        character(*), intent(in) :: name
        logical,      intent(in) :: val

        if (val) then
            call print_key_value_pair_character(name, 'true')
        else
            call print_key_value_pair_character(name, 'false')
        endif
    end subroutine print_key_value_pair_logical

    subroutine print_key_value_pair_character(name, val)
        character(*), intent(in) :: name
        character(*), intent(in) :: val
        character(40)            :: fix_length_name

        fix_length_name = name
        write (*, "(a, '=', a15)") fix_length_name, val
    end subroutine print_key_value_pair_character

end module iomanip
