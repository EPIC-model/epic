program free_slip
    use hdf5
    use constants, only : pi
    use parcel_container
    use interpolation, only : par2grid
    use options, only : parcel_info, grid
    use parameters, only : extent, lower, update_parameters, vcell
    use writer
    implicit none

    double precision :: volume_f(0:5, -1:6, 1)
    integer :: i, j
    integer(hid_t) :: group, step_group
    double precision, parameter :: angle = 0.5 * pi
    double precision, parameter :: lam = 2.0
    character(:), allocatable :: step
    character(:), allocatable :: name

    grid = (/5, 5/)
    extent =  (/0.5, 0.5/)
    lower = (/-0.25, -0.25/)

    call update_parameters()

    call parcel_alloc(25)

    do j = 1, 5
        do i = 1, 5
            parcels%position(i + 5 * (j-1), 1) = -0.3 + 0.1 * i
            parcels%position(i + 5 * (j-1), 2) = -0.2 + 0.1 * (j - 1)
        enddo
    enddo
    n_parcels = 25

    volume_f = 0.0

    parcel_info%is_elliptic = .true.

    parcels%volume = 0.05 * vcell

    ! b11
    parcels%B(:, 1) = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2

    ! b12
    parcels%B(:, 2) = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)



    call par2grid(parcels, parcels%volume, volume_f)

    call open_h5_file('free_slip.hdf5')

    call write_h5_double_scalar_step_attrib(0, "t", 0.0d0)

    call write_h5_double_scalar_step_attrib(0, "dt", 0.0d0)
    call write_h5_integer_scalar_step_attrib(0, "num parcel", n_parcels)

    call write_h5_parcels(0)

    !!!
    step = trim(get_step_group_name(0))
    name = step // "/fields"
    step_group = open_h5_group(step)
    group = open_h5_group(name)

    call write_h5_dataset_3d(name, "volume", &
                             volume_f(:, 0:, :))

    call h5gclose_f(group, h5err)
    call h5gclose_f(step_group, h5err)
    !!!


    group = open_h5_group("/")
    call write_h5_integer_scalar_attrib(group, "nsteps", 1)
    call h5gclose_f(group, h5err)

    group = open_h5_group("parcel")
    call write_h5_integer_scalar_attrib(group, "n_per_cell", 1)
    call write_h5_logical_attrib(group, "is_random", .false.)
    call write_h5_integer_scalar_attrib(group, "seed", 42)
    call write_h5_logical_attrib(group, "is_elliptic", .true.)
    call write_h5_double_scalar_attrib(group, "lambda", 3.0d0)
    call h5gclose_f(group, h5err)

    group = open_h5_group("mesh")
    call write_h5_double_vector_attrib(group, "extent", extent)
    call write_h5_double_vector_attrib(group, "origin", lower)
    call write_h5_integer_vector_attrib(group, "grid", grid)
    call h5gclose_f(group, h5err)


    call close_h5_file


    call parcel_dealloc()


end program free_slip
