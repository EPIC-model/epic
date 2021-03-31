module parameters
    use types, only : mesh_type
    implicit none
    !
    ! global parameters
    !

    ! print more info if true
    logical :: verbose = .false.

    ! configuration file
    character(len=32) :: filename = ''

    !
    ! output parameters
    !
    type output_info_type
        integer                     :: h5freq   = 1
        character(len=32)           :: h5fname  = ''
    end type output_info_type

    type(output_info_type) :: output

    !
    ! domain info
    !

    type(mesh_type) :: mesh

    !
    ! parcel parameters
    !
    type parcel_info_type
        integer          :: n_per_cell   = 4           ! number of parcels per cell (need to be a square)
        logical          :: is_random    = .false.     ! random parcel initialization
        integer          :: seed         = 42          ! seed of random initialization
        logical          :: is_elliptic  = .true.      ! use elliptic model
        double precision :: lambda       = 5.0         ! max. ellipse aspect ratio a/b
        integer          :: split_freq   = 1           ! split frequency, 1 = call split subroutine every step
        character(10)    :: merge_type   = 'geometric' ! merge method in use (other option: 'optimal')
        integer          :: merge_freq   = 1           ! merge frequency, 1 = call merge subroutine every step
        double precision :: vmin         = 0.001       ! min. allowed parcel volume
    end type parcel_info_type

    type(parcel_info_type) :: parcel_info

    !
    ! stepper parameters
    !

    ! time limit
    type time_info_type
        double precision :: limit       = 0.0       ! time limit
        double precision :: dt          = 0.0       ! time step
        logical          :: is_adaptive = .false.
    end type time_info_type

    type(time_info_type) :: time

    !
    ! interpolation
    !
    character(len=32) :: interpl = 'trilinear'

    !
    ! TaylorGreen flow parameters
    !
    type flow_type
        double precision :: amp(2) = 1.0    ! amplitudes
        double precision :: freq(2) = 1.0   ! frequencies
        double precision :: phase(2) = 1.0  ! phase shift
    end type flow_type

    type(flow_type) :: flow

end module parameters
