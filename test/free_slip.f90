program free_slip
    use hdf5
    use constants, only : pi
    use parcel_container
    use writer
    implicit none

    double precision :: volume_f(5, 3)
    integer :: i, j
    integer(hid_t) :: group
    double precision :: origin(2) = (/-0.3, -0.25/)
    double precision :: extent(2) =  (/0.6, 0.5/)
    integer :: grid(2) = (/5, 2/)
    double precision, parameter :: angle = 0.5 * pi
    double precision, parameter :: lam = 3.0
    double precision :: B22, evec(2)

    call parcel_alloc(10)

    do j = 1, 2
        do i = 1, 5
            parcels%position(i + 5 * (j-1), 1) = -0.25 + 0.1 * i
            parcels%position(i + 5 * (j-1), 2) = 0.2 * (j - 1)
        enddo
    enddo
    n_parcels = 10



    parcels%volume = 0.25 * product(extent / (grid - 1))

    ! b11
    parcels%B(:, 1) = lam * cos(angle) ** 2 + 1.0 / lam * sin(angle) ** 2

    ! b12
    parcels%B(:, 2) = 0.5 * (lam - 1.0 / lam) * sin(2.0 * angle)

    print *, angle
    B22 = lam * sin(angle) ** 2 + 1.0 / lam * cos(angle) ** 2

    evec(1) = lam - B22
    evec(2) = parcels%B(1, 2)

    print *, "before:", evec

    if (abs(evec(1)) <= epsilon(lam) .and. evec(1) /= 0.0d0) then
        evec(1) = evec(1) + epsilon(evec(1))
    else
        print *, "else"
        evec(2) = evec(2) + epsilon(evec(2))
    endif


    evec = evec / norm2(evec)
    print *, "here:", evec, norm2(evec)



    call open_h5_file('free_slip.hdf5')

    call write_h5_double_scalar_step_attrib(0, "t", 0.0d0)

    call write_h5_double_scalar_step_attrib(0, "dt", 0.0d0)
    call write_h5_integer_scalar_step_attrib(iter, "num parcel", n_parcels)

    call write_h5_parcels(0)


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
    call write_h5_double_vector_attrib(group, "origin", origin)
    call write_h5_integer_vector_attrib(group, "grid", grid)
    call h5gclose_f(group, h5err)


    call close_h5_file





    call parcel_dealloc()


end program free_slip
