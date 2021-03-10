module options
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.

    ! configuration file
    character(len=32) :: filename = ''

    !
    ! output options
    !
    integer :: h5freq = 1

    !
    ! parcel options
    !

    ! maximum number of allowed parcels
    integer :: max_num_parcels = 1e4

    ! number of parcels per cell
    integer :: n_per_cell = 4

    ! random parcel initialization
    logical :: is_random = .false.

    ! seed of random initialization
    integer :: seed = 42

    ! use elliptic model
    logical :: is_elliptic = .true.

    !
    ! stepper options
    !

    ! time limit
    double precision :: tmax = 0.0

    ! time step
    double precision :: dt = 0.0

    ! adaptive time step
    logical :: is_adaptive = .false.

end module options
