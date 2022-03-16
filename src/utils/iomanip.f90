module iomanip

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

end module iomanip
