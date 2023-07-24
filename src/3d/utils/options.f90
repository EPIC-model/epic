! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    use constants, only : zero, one, two, pi, four, twopi
    use netcdf_writer
    use mpi_utils, only : mpi_stop
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.

    ! if a restarted simulation
    logical :: l_restart = .false.

    ! configuration file
    character(len=512) :: filename = ''

    ! restart file
    character(len=512) :: restart_file = ''

    ! field input file
    character(len=512)  :: field_file = ''

    ! ls rk order
    integer :: rk_order = 4


    type bndry_info
        ! the rms of the boundary zeta must be smaller than
        ! this threshold times the rms of the interior
        ! zeta in order to keep zeta = 0 and dzeta/dt = 0
        ! on the boundary;
        ! each boundary is checked independently
        double precision :: zeta_tol = 1.0e-3

        ! this option is enabled, it makes the code to ignore the
        ! zeta boundary flag; hence a simulation can develop a
        ! non-zero vertical vorticity component (zeta) on both
        ! boundaries. A warning is printed to the standard
        ! output if it is enabled.
        ! WARNING: You need to be aware of the consequences
        !          on the physics when enabling this option!
        logical :: l_ignore_bndry_zeta_flag = .false.
    end type bndry_info

    type(bndry_info) :: boundary

    !
    ! output options
    !
    type info
        double precision                  :: field_freq         = one
        logical                           :: write_fields       = .true.
        character(len=64), dimension(128) :: field_list         = ''
        double precision                  :: parcel_freq        = one
        logical                           :: overwrite          = .false.
        logical                           :: write_parcels      = .true.
        double precision                  :: parcel_stats_freq  = one
        logical                           :: write_parcel_stats = .true.
        double precision                  :: field_stats_freq   = one
        logical                           :: write_field_stats  = .true.
        character(len=512)                :: basename           = ''
    end type info

    type(info) :: output

    !
    ! domain options
    !
    logical :: allow_larger_anisotropy = .false.

    !
    ! parcel options
    !
    type parcel_type
        double precision :: size_factor      = 1.0d0    ! factor to increase max. number of parcels
        double precision :: grow_factor      = 1.2d0    ! factor to increase the parcel container size
                                                        ! in the parcel splitting routine
        double precision :: shrink_factor    = 0.8d0    ! factor to reduce the parcel container size
        integer          :: n_per_cell       = 8        ! number of parcels per cell (need to be a cube)
        double precision :: lambda_max       = four     ! max. ellipse aspect ratio a/b
        double precision :: min_vratio       = 40.0d0   ! minimum ratio of grid cell volume / parcel volume
        integer          :: correction_iters = 2        ! parcel correction iterations
        double precision :: gradient_pref    = 1.8d0    ! prefactor for gradient descent
        double precision :: max_compression  = 0.5d0    ! parameter for gradient descent
                                                        ! (limits the shift in parcel position)
    end type parcel_type

    type(parcel_type) :: parcel

    ! time limit
    type time_info_type
        double precision :: initial     = zero       ! initial time
        double precision :: limit       = zero       ! time limit
        double precision :: alpha       = 0.2d0      ! factor for adaptive time stepping with strain and buoyancy
                                                     ! gradient
        logical          :: precise_stop = .false.   ! stop at the exact limit
    end type time_info_type

    type(time_info_type) :: time

    contains
        ! parse configuration file
        ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
        subroutine read_config_file
            integer :: ios
            integer :: fn = 1
            logical :: exists = .false.

            ! namelist definitions
            namelist /EPIC/ field_file, rk_order, boundary, output, parcel, time

            ! check whether file exists
            inquire(file=filename, exist=exists)

            if (exists .eqv. .false.) then
                call mpi_stop(&
                    'Error: input file "' // trim(filename) // '" does not exist.')
            endif

            ! open and read Namelist file.
            open(action='read', file=filename, iostat=ios, newunit=fn)

            read(nml=EPIC, iostat=ios, unit=fn)

            if (ios /= 0) then
                call mpi_stop('Error: invalid Namelist format.')
            end if

            close(fn)

            ! check whether NetCDF files already exist
            inquire(file=output%basename, exist=exists)

            if (exists) then
                call mpi_stop(&
                    'Error: output file "' // trim(output%basename) // '" already exists.')
            endif

        end subroutine read_config_file

        subroutine write_netcdf_options(ncid)
            integer, intent(in) :: ncid

#ifdef ENABLE_VERBOSE
            call write_netcdf_attribute(ncid, "verbose", verbose)
#endif
            call write_netcdf_attribute(ncid, "field_file", field_file)

            call write_netcdf_attribute(ncid, "rk_order", rk_order)

            call write_netcdf_attribute(ncid, "zeta_tol", boundary%zeta_tol)
            call write_netcdf_attribute(ncid, "l_ignore_bndry_zeta_flag", &
                                        boundary%l_ignore_bndry_zeta_flag)

            call write_netcdf_attribute(ncid, "allow_larger_anisotropy", &
                                               allow_larger_anisotropy)

            call write_netcdf_attribute(ncid, "size_factor", parcel%size_factor)
            call write_netcdf_attribute(ncid, "n_per_cell", parcel%n_per_cell)
            call write_netcdf_attribute(ncid, "lambda_max", parcel%lambda_max)
            call write_netcdf_attribute(ncid, "min_vratio", parcel%min_vratio)
            call write_netcdf_attribute(ncid, "correction_iters", parcel%correction_iters)
            call write_netcdf_attribute(ncid, "gradient_pref", parcel%gradient_pref)
            call write_netcdf_attribute(ncid, "max_compression", parcel%max_compression)

            call write_netcdf_attribute(ncid, "parcel_freq", output%parcel_freq)
            call write_netcdf_attribute(ncid, "field_freq", output%field_freq)
            call write_netcdf_attribute(ncid, "parcel_stats_freq", output%parcel_stats_freq)
            call write_netcdf_attribute(ncid, "write_parcel_stats", output%write_parcel_stats)
            call write_netcdf_attribute(ncid, "field_stats_freq", output%field_stats_freq)
            call write_netcdf_attribute(ncid, "write_field_stats", output%write_field_stats)
            call write_netcdf_attribute(ncid, "write_fields", output%write_fields)
            call write_netcdf_attribute(ncid, "overwrite", output%overwrite)
            call write_netcdf_attribute(ncid, "write_parcels", output%write_parcels)
            call write_netcdf_attribute(ncid, "basename", trim(output%basename))

            call write_netcdf_attribute(ncid, "limit", time%limit)
            call write_netcdf_attribute(ncid, "initial", time%initial)
            call write_netcdf_attribute(ncid, "precise_stop", time%precise_stop)
            call write_netcdf_attribute(ncid, "alpha", time%alpha)

        end subroutine write_netcdf_options

end module options
