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
    character(len=512) :: filename = ''

    ! field input file
    character(len=512)  :: field_file = ''
    double precision    :: field_tol  = 1.0d-10

    !
    ! output options
    !
    type h5_info
        integer             :: h5_field_freq         = 1
        logical             :: h5_write_fields       = .true.
        integer             :: h5_parcel_freq        = 1
        logical             :: h5_overwrite          = .false.
        logical             :: h5_write_parcels      = .true.
        integer             :: h5_diagnostics_freq   = 1
        logical             :: h5_write_diagnostics  = .true.
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
        integer          :: n_per_cell   = 4              ! number of parcels per cell (need to be a square)
        double precision :: lambda_max   = five           ! max. ellipse aspect ratio a/b
        integer          :: split_freq   = 1              ! split frequency, 1 = call split subroutine every step
        character(15)    :: merge_type   = 'bi-geometric' ! merge method in use (other option: 'optimal')
        integer          :: merge_freq   = 1              ! merge frequency, 1 = call merge subroutine every step
        double precision :: vmin_fraction   = 36.0d0      ! minimum parcel area fraction (in merge)
        integer          :: correction_freq = 1           ! parcel correction frequency, 1 = call module every step
        integer          :: correction_iters= 1           ! parcel correction iterations
        logical          :: apply_laplace = .true.        ! use Laplacian parcel correction
        logical          :: apply_gradient = .true.       ! use gradient descent to adjust parcel positions on small scale
        double precision :: gradient_pref= 1.8d0          ! prefactor for gradient descent
        double precision :: max_compression= 0.5d0        ! parameter for gradient descent (limits the shift in parcel position)
        double precision :: vmax_fraction = 2.0           ! maximum parcel area fraction

    end type parcel_type

    type(parcel_type) :: parcel

    ! time limit
    type time_info_type
        double precision :: limit       = zero       ! time limit
        double precision :: dt          = zero       ! time step
        logical          :: is_adaptive = .false.
        double precision :: alpha_s     = 0.1d0      ! factor for adaptive time stepping with strain
        double precision :: alpha_b     = 0.1d0      ! factor for adaptive time stepping with buoyancy gradient
    end type time_info_type

    type(time_info_type) :: time

end module options
