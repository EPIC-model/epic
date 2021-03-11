module fields
    use parameters, only : mesh
    implicit none

    double precision, allocatable, dimension(:, :, :) :: &
        velocity

    contains
        function get_mesh_spacing() result(dx)
            double precision :: dx(2)
            dx = mesh%extent / (mesh%grid - 1)
        end function get_mesh_spacing

end module fields
