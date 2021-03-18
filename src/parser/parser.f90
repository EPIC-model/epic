module parser
    use parameters
    use hdf5
    use writer
    implicit none

    private
    public :: read_config_file, &
              write_h5_params

    type(parcel_info_type) :: parcel

    contains

        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /MODEL/ output, mesh, parcel, time, flow, interpl

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                print *, 'Error: input file "', trim(filename), '" does not exist.'
                stop
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=MODEL, iostat=ios, unit=fn)

            if (ios /= 0) then
                print *, 'Error: invalid Namelist format.'
                stop
            end if

            close(fn)

            ! check whether h5 file already exists
            inquire(file=output%h5fname, exist=exists)

            if (exists) then
                print *, 'Error: output file "', trim(output%h5fname), '" already exists.'
                stop
            endif

            ! update parcel parameters
            parcel_info = parcel

        end subroutine read_config_file

        subroutine write_h5_params
            integer(hid_t) :: group

            call open_h5_file(trim(output%h5fname))

            !
            ! write parcel info
            !
            group = open_h5_group("parcel")

            call write_h5_integer_scalar_attrib(group, "n_per_cell", parcel_info%n_per_cell)
            call write_h5_logical_attrib(group, "is_random", parcel_info%is_random)
            call write_h5_integer_scalar_attrib(group, "seed", parcel_info%seed)
            call write_h5_logical_attrib(group, "is_elliptic", parcel_info%is_elliptic)

            call h5gclose_f(group, h5err)

            !
            ! write output info
            !
            group = open_h5_group("output")

            call write_h5_integer_scalar_attrib(group, "h5freq", output%h5freq)

            call h5gclose_f(group, h5err)

            !
            ! write stepper info
            !
            group = open_h5_group("time")

            call write_h5_double_scalar_attrib(group, "limit", time%limit)
            call write_h5_logical_attrib(group, "is_adaptive", time%is_adaptive)

            call h5gclose_f(group, h5err)

            !
            ! mesh info
            !
            group = open_h5_group("mesh")

            call write_h5_double_vector_attrib(group, "extent", mesh%extent)
            call write_h5_double_vector_attrib(group, "origin", mesh%origin)
            call write_h5_integer_vector_attrib(group, "grid", mesh%grid)
            call write_h5_character_vector_attrib(group, "bc", mesh%bc)

            call h5gclose_f(group, h5err)

            call close_h5_file

        end subroutine write_h5_params

end module parser
