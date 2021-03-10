module parser
    use options
    implicit none

    private
    public :: read_config_file

    type :: mesh_info_type
        double precision  :: origin(2)
        double precision  :: extent(2)
        integer           :: grid(2)
        character(len=16) :: bc(2)
    end type mesh_info_type

    type(mesh_info_type) :: mesh

    type parcel_info_type
        integer :: n_per_cell
        logical :: is_random
        integer :: seed
        logical :: is_elliptic
    end type parcel_info_type

    type(parcel_info_type) :: parcel

    type time_info_type
        double precision :: tmax        ! time limit
        double precision :: dt          ! time step
        logical          :: is_adaptive
    end type time_info_type


    type(time_info_type) :: time

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

            call update_options

        end subroutine read_config_file

        ! after parsing we need to update the global options
        subroutine update_options
            ! update parcel options
            n_per_cell  = parcel%n_per_cell
            is_random   = parcel%is_random
            seed        = parcel%seed
            is_elliptic = parcel%is_elliptic

            ! update stepper options
            tmax        = time%tmax
            dt          = time%dt
            is_adaptive = time%is_adaptive
        end subroutine update_options


end module parser
