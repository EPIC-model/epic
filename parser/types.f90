module types
    implicit none

    type :: mesh_type
        double precision  :: origin(2)
        double precision  :: extent(2)
        integer           :: grid(2)
        character(len=16) :: bc(2)
    end type mesh_type

    type time_info_type
        double precision :: limit       ! time limit
        double precision :: dt          ! time step
        logical          :: adaptive
    end type time_info_type

    type parcel_info_type
        integer :: n_per_cell
        logical :: is_random
        integer :: seed
        logical :: is_elliptic
    end type parcel_info_type

end module types
