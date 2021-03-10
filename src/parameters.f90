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

    ! random parcel initialization
    logical :: is_random = .false.

    ! seed of random initialization
    integer :: seed = 42

    ! use elliptic model
    logical :: is_elliptic = .true.

    !
    ! stepper parameters
    !

    ! time limit
    double precision :: tmax = 0.0

    ! time step
    double precision :: dt = 0.0

    ! adaptive time step
    logical :: is_adaptive = .false.

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
