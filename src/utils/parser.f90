! =============================================================================
!             This module to parses the user specifications provided
!             by the namelist.
! =============================================================================
module parser
    use constants
    use options
    implicit none

    contains

        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /EPIC/ field_file, field_tol, output, parcel, time

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                print *, 'Error: input file "', trim(filename), '" does not exist.'
                stop
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=EPIC, iostat=ios, unit=fn)

            if (ios /= 0) then
                print *, 'Error: invalid Namelist format.'
                stop
            end if

            close(fn)

            ! check whether h5 files already exist
            inquire(file=output%h5_basename, exist=exists)

            if (exists) then
                print *, 'Error: output file "', trim(output%h5_basename), '" already exists.'
                stop
            endif

        end subroutine read_config_file
end module parser
