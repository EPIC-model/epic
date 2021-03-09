module types
    implicit none

    type :: mesh_type
        double precision  :: origin(2)
        double precision  :: extent(2)
        integer           :: grid(2)
        character(len=16) :: bc(2)
    end type mesh_type


end module types
