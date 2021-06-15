! =============================================================================
!             This module to parses the user specifications provided
!             by the namelist.
! =============================================================================
module parser
    use constants
    use options
    use h5_writer
    implicit none

    contains

        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /EPIC/ model, output, parcel, stepper, time

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


        subroutine write_h5_options(h5fname)
            character(*), intent(in) :: h5fname
            integer(hid_t)           :: h5handle, gopts, group, gexec
            character(len=8)         :: date, tmp2
            character(len=10)        :: ctime, tmp1

            call open_h5_file(h5fname, H5F_ACC_RDWR_F, h5handle)

            call create_h5_group(h5handle, "execution", gexec)
                ! 15 June 2021
                ! https://gcc.gnu.org/onlinedocs/gfortran/DATE_005fAND_005fTIME.html
                call date_and_time(date=date, time=ctime)
                ! 15 June 2021
                ! https://stackoverflow.com/questions/13755762/access-character-at-specific-index-in-a-string-in-fortran
                tmp1 = date(1:4) // '-' // date(5:6) // '-' // date(7:8)
                tmp2 = ctime(1:2) // '-' // ctime(3:4) // '-' // ctime(5:6)
                call write_h5_char_scalar_attrib(gexec, "date", tmp1)
                call write_h5_char_scalar_attrib(gexec, "time", tmp2)
            call close_h5_group(gexec)

            call create_h5_group(h5handle, "options", gopts)

            call create_h5_group(gopts, "parcel", group)
                call write_h5_int_scalar_attrib(group, "n_per_cell", parcel%n_per_cell)
                call write_h5_logical_attrib(group, "is_random", parcel%is_random)
                call write_h5_int_scalar_attrib(group, "seed", parcel%seed)
                call write_h5_logical_attrib(group, "is_elliptic", parcel%is_elliptic)
                call write_h5_double_scalar_attrib(group, "lambda", parcel%lambda)
                call write_h5_int_scalar_attrib(group, "h5_parcel_freq", output%h5_parcel_freq)
            call close_h5_group(group)

            call create_h5_group(gopts, "time", group)
                call write_h5_double_scalar_attrib(group, "time limit", time%limit)
                call write_h5_logical_attrib(group, "is_adaptive", time%is_adaptive)
            call close_h5_group(group)

            call create_h5_group(gopts, "stepper", group)
                call write_h5_char_scalar_attrib(group, "method", stepper)
            call close_h5_group(group)

            call close_h5_group(gopts)
            call close_h5_file(h5handle)
        end subroutine write_h5_options

end module parser
