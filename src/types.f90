module types
    implicit none

    type :: mesh_type
        double precision  :: origin(2)  ! mesh origin
        double precision  :: extent(2)  ! mesh extent
        integer           :: grid(2)    ! number of grid points per dimension
    end type mesh_type
end module types
