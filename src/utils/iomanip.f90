module iomanip

    interface print_quantity
        module procedure :: print_quantity_double
        module procedure :: print_quantity_integer
        module procedure :: print_quantity_logical
        module procedure :: print_quantity_character
    end interface print_quantity

    private :: print_quantity_double,     &
               print_quantity_integer,    &
               print_quantity_logical,    &
               print_quantity_character
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


    subroutine print_quantity_double(name, val, unit)
        character(*),           intent(in) :: name
        double precision,       intent(in) :: val
        character(*), optional, intent(in) :: unit
        character(64)                      :: fix_length_name
        character(16)                      :: fix_length_unit = ''

        fix_length_name = name
        if (present(unit)) then
            fix_length_unit = unit
        endif
        write (*, "(a, 1p,e14.7, a)") fix_length_name, val, fix_length_unit
    end subroutine print_quantity_double

    subroutine print_quantity_integer(name, val, unit)
        character(*),           intent(in) :: name
        integer,                intent(in) :: val
        character(*), optional, intent(in) :: unit
        character(64)                      :: fix_length_name
        character(16)                      :: fix_length_unit = ''

        fix_length_name = name
        if (present(unit)) then
            fix_length_unit = unit
        endif
        write (*, "(a, I14, a)") fix_length_name, val, fix_length_unit
    end subroutine print_quantity_integer

    subroutine print_quantity_logical(name, val, unit)
        character(*),           intent(in) :: name
        logical,                intent(in) :: val
        character(*), optional, intent(in) :: unit

        if (val) then
            call print_quantity_character(name, 'true')
        else
            call print_quantity_character(name, 'false')
        endif
    end subroutine print_quantity_logical

    subroutine print_quantity_character(name, val, unit)
        character(*),           intent(in) :: name
        character(*),           intent(in) :: val
        character(*), optional, intent(in) :: unit
        character(64)                      :: fix_length_name
        character(16)                      :: fix_length_unit = ''

        fix_length_name = name
        write (*, "(a, a14, a)") fix_length_name, val, fix_length_unit
    end subroutine print_quantity_character

end module iomanip
