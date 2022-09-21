module netcdf_dimensions
    use config, only : package
    implicit none

    character(1), parameter :: x_dim_name = 'x'
    character(1), parameter :: y_dim_name = 'y'
    character(1), parameter :: z_dim_name = 'z'
    character(1), parameter :: t_dim_name = 't'

    character(1), parameter :: x_axis_name = 'X'
    character(1), parameter :: y_axis_name = 'Y'
    character(1), parameter :: z_axis_name = 'Z'
    character(1), parameter :: t_axis_name = 'T'

    character(*), parameter :: version_name = package//'_version'

end module netcdf_dimensions
