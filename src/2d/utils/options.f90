! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    use constants, only : zero, one, two, pi, four
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.

    ! configuration file
    character(len=512) :: filename = ''

    ! field input file
    character(len=512)  :: field_file = ''
    double precision    :: field_tol  = 1.0d-10

    !
    ! output options
    !
    type h5_info
        double precision    :: h5_field_freq         = one
        logical             :: h5_write_fields       = .true.
        double precision    :: h5_parcel_freq        = one
        logical             :: h5_overwrite          = .false.
        logical             :: h5_write_parcels      = .true.
        double precision    :: h5_parcel_stats_freq  = one
        logical             :: h5_write_parcel_stats = .true.
        double precision    :: h5_field_stats_freq   = one
        logical             :: h5_write_field_stats  = .true.
        character(len=512)  :: h5_basename           = ''
    end type h5_info

    type(h5_info) :: output

    !
    ! domain options
    !
    logical :: allow_larger_anisotropy = .false.

    !
    ! parcel options
    !
    type parcel_type
        integer          :: n_per_cell       = 9        ! number of parcels per cell (need to be a square)
        double precision :: lambda_max       = four     ! max. ellipse aspect ratio a/b
        double precision :: min_vratio       = 40.0d0   ! minimum ratio of grid cell volume / parcel volume
        integer          :: correction_iters = 2        ! parcel correction iterations
        double precision :: gradient_pref    = 1.8d0    ! prefactor for gradient descent
        double precision :: max_compression  = 0.5d0    ! parameter for gradient descent (limits the shift in parcel position)
        double precision :: max_vratio       = 2.89     ! maximum ratio of grid cell volume / parcel volume

    end type parcel_type

    type(parcel_type) :: parcel

    ! time limit
    type time_info_type
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

end module options
