! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    use constants, only : zero, one, two, pi, five
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.

    ! configuration file
    character(len=32) :: filename = ''

    character(len=32) :: model = 'TaylorGreen'

    !
    ! output options
    !
    type output_info_type
        integer                     :: h5freq   = 1
        character(len=32)           :: h5fname  = ''
    end type output_info_type

    type(output_info_type) :: output

    !
    ! domain options
    !
    type box_type
        integer          :: nc(2) = (/20, 20/)                  ! number of cells
        double precision :: extent(2) = (/pi, pi/)
    end type box_type

    type(box_type) :: box

    !
    ! parcel options
    !
    type parcel_info_type
        integer          :: n_per_cell   = 4              ! number of parcels per cell (need to be a square)
        logical          :: is_random    = .false.        ! random parcel initialization
        integer          :: seed         = 42             ! seed of random initialization
        logical          :: is_elliptic  = .true.         ! use elliptic model
        double precision :: lambda       = five           ! max. ellipse aspect ratio a/b
        integer          :: split_freq   = 1              ! split frequency, 1 = call split subroutine every step
        character(15)    :: merge_type   = 'bi-geometric' ! merge method in use (other option: 'optimal')
        integer          :: merge_freq   = 1              ! merge frequency, 1 = call merge subroutine every step
        double precision :: vfraction    = 36.0d0         ! volume fraction in merge
        integer          :: correction_freq = 1           ! parcel correction frequency, 1 = call module every step
        integer          :: correction_iters= 1           ! parcel correction iterations
        logical          :: apply_laplace = .true.        ! use Laplacian parcel correction
        logical          :: apply_gradient = .true.       ! use gradient descent to adjust parcel positions on small scale
        double precision :: gradient_pref= 0.3d0          ! prefactor for gradient descent
        double precision :: vmaxfraction = 2.0         ! prefactor for gradient descent

    end type parcel_info_type

    type(parcel_info_type) :: parcel_info

    !
    ! stepper options
    !
    character(len=16) :: stepper = 'classic-rk4'

    ! time limit
    type time_info_type
        double precision :: limit       = zero       ! time limit
        double precision :: dt          = zero       ! time step
        logical          :: is_adaptive = .false.
        double precision :: alpha       = 0.025d0   ! factor for adaptive time stepping with strain
        double precision :: dt_max      = 0.125d0   ! maximum time step
    end type time_info_type

    type(time_info_type) :: time

    !
    ! interpolation
    !
    character(len=32) :: interpl = 'trilinear'

    !
    ! TaylorGreen flow options
    !
    type flow_type
        double precision :: amp(2) = one    ! amplitudes
        double precision :: freq(2) = one   ! frequencies
        double precision :: phase(2) = one  ! phase shift
    end type flow_type

    type(flow_type) :: flow

end module options
