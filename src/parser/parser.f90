module parser
    use parameters
    implicit none

    private
    public :: read_config_file

    type(parcel_info_type) :: parcel

    contains

        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer                                :: ios
            integer                                :: fn = 1

            ! namelist definitions
            namelist /MODEL/ h5freq, mesh, parcel, time

            ! check whether file exists.
            inquire(file=filename, iostat=ios)

            if (ios /= 0) then
                write (*, *) 'Error: input file "', filename, '" does not exist.'
                return
            endif

            ! open and read Namelist file.
            open (action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=MODEL, iostat=ios, unit=fn)

            if (ios /= 0) then
                write (*, *) 'Error: invalid Namelist format.'
            end if

            close(fn)

            ! update parcel parameters
            parcel_info = parcel

        end subroutine read_config_file

end module parser
