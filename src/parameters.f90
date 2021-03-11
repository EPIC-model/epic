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
    integer :: h5freq = 1

    !
    ! domain info
    !

    type(mesh_type) :: mesh

    !
    ! parcel parameters
    !
    type parcel_info_type
        logical :: is_random    = .false.   ! random parcel initialization
        integer :: seed         = 42        ! seed of random initialization
        logical :: is_elliptic  = .true.    ! use elliptic model
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
    ! TaylorGreen flow parameters
    !
    ! amplitudes
    double precision :: amp(2) = 1.0

    ! frequencies
    double precision :: freq(2) = 1.0

    ! phase shift
    double precision :: phase(2) = 1.0

end module parameters
