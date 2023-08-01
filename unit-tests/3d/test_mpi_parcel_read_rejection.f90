! =============================================================================
!          Test parallel parcel reading using the rejection method.
!
!   If a parcel file does not contain the start_index dataset,
!   all MPI ranks read all parcels, but reject parcels that do not
!   belong to the sub-domain owned by the rank.
! =============================================================================
program test_mpi_parcel_read_rejection
    use unit_test
    use options, only : parcel
    use constants, only : zero, f12
    use parcel_container
    use parcel_netcdf
    use mpi_environment
    use mpi_layout
    use parameters, only : lower, update_parameters, extent, nx, ny, nz, dx, max_num_parcels
    use mpi_timer
    implicit none

    logical            :: passed = .true.
    integer, parameter :: n_per_dim = 2
    double precision   :: res
    integer            :: i, j, k, l, ix, iy, iz, n_parcels_before
    double precision   :: im, corner(3)
    double precision   :: x_sum, y_sum, z_sum

    integer            :: ncid
    integer            :: npar_dim_id, vol_id, theta_id,      &
                          x_pos_id, y_pos_id, z_pos_id,     &
                          x_vor_id, y_vor_id, z_vor_id,     &
                          b11_id, b12_id, b13_id,           &
                          b22_id, b23_id,                   &
                          t_axis_id, t_dim_id, mpi_dim_id
#ifndef ENABLE_DRY_MODE
    integer :: qv_id, ql_id
#endif

    character(len=512) :: ncbasename

    character(len=512) :: ncfname

    call mpi_env_initialise

    passed = (passed .and. (world%err == 0))

    nx = 32
    ny = 32
    nz = 32
    lower  = zero
    extent =  one

    ! needed for max_num_parcels
    parcel%min_vratio = 40.d0
    parcel%size_factor = 1

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call register_timer('parcel I/O', parcel_io_timer)

    !
    ! write parcels first
    !

    n_parcels = n_per_dim ** 3 * nz * (box%hi(2)-box%lo(2)+1) * (box%hi(1)-box%lo(1)+1)
    n_parcels_before = n_parcels
    n_total_parcels = n_per_dim ** 3 * nz * ny * nx
    call parcel_alloc(max_num_parcels)

    im = one / dble(n_per_dim)

    l = 1
    do iz = 0, nz-1
        do iy = box%lo(2), box%hi(2)
            do ix = box%lo(1), box%hi(1)
                corner = lower + dble((/ix, iy, iz/)) * dx
                do k = 1, n_per_dim
                    do j = 1, n_per_dim
                        do i = 1, n_per_dim
                            parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                            parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                            parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                            l = l + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    ! we cannot check individual positions since the order of the parcels can be arbitrary
    x_sum = sum(parcels%position(1, 1:n_parcels))
    y_sum = sum(parcels%position(2, 1:n_parcels))
    z_sum = sum(parcels%position(3, 1:n_parcels))

    parcels%B(:, 1:n_parcels) = world%rank + 1
    parcels%volume(1:n_parcels) = world%rank + 1
    parcels%vorticity(:, 1:n_parcels) = world%rank + 1
    parcels%theta(1:n_parcels) = world%rank + 1

#ifndef ENABLE_DRY_MODE
    parcels%qv(1:n_parcels) = world%rank + 1
    parcels%ql(1:n_parcels) = world%rank + 1
#endif

    call create_file('nctest')

    passed = (passed .and. (ncerr == 0))

    call write_parcels(t = 10.0d0)

    passed = (passed .and. (ncerr == 0))

    !
    ! read parcels now
    !

    parcels%position = 0
    parcels%B = 0
    parcels%volume = 0
    parcels%vorticity = 0
    parcels%theta = 0
#ifndef ENABLE_DRY_MODE
    parcels%qv = 0
    parcels%ql = 0
#endif

    call read_netcdf_parcels('nctest_0000000001_parcels.nc')

    passed = (passed .and. (n_parcels == n_parcels_before))

    if (passed) then
        res = dble(world%rank + 1)

        passed = (passed .and. (abs(sum(parcels%position(1, 1:n_parcels)) - x_sum) == zero))
        passed = (passed .and. (abs(sum(parcels%position(2, 1:n_parcels)) - y_sum) == zero))
        passed = (passed .and. (abs(sum(parcels%position(3, 1:n_parcels)) - z_sum) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(1, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(2, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(3, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(4, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(5, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%volume(1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(1, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(2, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(3, 1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%theta(1:n_parcels) - res)) == zero))
#ifndef ENABLE_DRY_MODE
        passed = (passed .and. (maxval(abs(parcels%qv(1:n_parcels) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%ql(1:n_parcels) - res)) == zero))
#endif
    endif

    call parcel_dealloc

    call delete_netcdf_file(ncfname='nctest_0000000001_parcels.nc')

    passed = (passed .and. (ncerr == 0))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI parcel read', passed)
    endif

    call mpi_env_finalise


    contains
        subroutine create_file(basename)
            character(*), intent(in)  :: basename
            integer                   :: dimids(2)

            ncfname =  basename // '_' // zfill(1) // '_parcels.nc'

            ncbasename = basename

            call create_netcdf_file(ncfname, .true., ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                       &
                                   version_tag=package_version,     &
                                   file_type='parcels',             &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_parcels',                  &
                                         dimsize=n_total_parcels,           &
                                         dimid=npar_dim_id)

            call define_netcdf_dimension(ncid=ncid,                 &
                                         name='world%size',           &
                                         dimsize=world%size,          &
                                         dimid=mpi_dim_id)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            dimids = (/npar_dim_id, t_dim_id/)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='x_position',                   &
                                       long_name='x position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=x_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='y_position',                   &
                                       long_name='y position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=y_pos_id)

            call define_netcdf_dataset(ncid=ncid,                           &
                                       name='z_position',                   &
                                       long_name='z position component',    &
                                       std_name='',                         &
                                       unit='m',                            &
                                       dtype=NF90_DOUBLE,                   &
                                       dimids=dimids,                       &
                                       varid=z_pos_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B11',                              &
                                       long_name='B11 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b11_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B12',                              &
                                       long_name='B12 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b12_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B13',                              &
                                       long_name='B13 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b13_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B22',                              &
                                       long_name='B22 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b22_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='B23',                              &
                                       long_name='B23 element of shape matrix', &
                                       std_name='',                             &
                                       unit='m^2',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=b23_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='volume',                           &
                                       long_name='parcel volume',               &
                                       std_name='',                             &
                                       unit='m^3',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=vol_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='x_vorticity',                      &
                                       long_name='x vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=x_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='y_vorticity',                      &
                                       long_name='y vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=y_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='z_vorticity',                      &
                                       long_name='z vorticity component',       &
                                       std_name='',                             &
                                       unit='1/s',                              &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=z_vor_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='theta',                            &
                                       long_name='parcel potential temperature',&
                                       std_name='',                             &
                                       unit='K',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=theta_id)

#ifndef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='qv',                               &
                                       long_name='parcel water vapour spec. hum.',&
                                       std_name='',                             &
                                       unit='kg/kg',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=qv_id)

            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='ql',                               &
                                       long_name='parcel liquid water spec. hum.',&
                                       std_name='',                             &
                                       unit='kg/kg',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=ql_id)
#endif

            call close_definition(ncid)

            call close_netcdf_file(ncid)

        end subroutine create_file


        subroutine write_parcels(t)
            double precision, intent(in) :: t
            integer                      :: cnt(2), start(2)
            integer                      :: recvcounts(world%size)
            integer                      :: sendbuf(world%size), start_index

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! after this operation all MPI ranks know their starting index
            recvcounts = 1
            start_index = 0
            sendbuf = 0
            sendbuf(world%rank+1:world%size) = n_parcels
            sendbuf(world%rank+1) = 0

            call MPI_Reduce_scatter(sendbuf, start_index, recvcounts, MPI_INT, MPI_SUM, world%comm, world%err)

            ! we need to increase the start_index by 1
            ! since the starting index in Fortran is 1 and not 0.
            start_index = start_index + 1

            start = (/ start_index, 1 /)
            cnt   = (/ n_parcels,   1 /)

            call write_netcdf_dataset(ncid, x_pos_id, parcels%position(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, y_pos_id, parcels%position(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, z_pos_id, parcels%position(3, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, b11_id, parcels%B(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b12_id, parcels%B(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b13_id, parcels%B(3, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b22_id, parcels%B(4, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, b23_id, parcels%B(5, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, vol_id, parcels%volume(1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, parcels%vorticity(1, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, parcels%vorticity(2, 1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, parcels%vorticity(3, 1:n_parcels), start, cnt)

            call write_netcdf_dataset(ncid, theta_id, parcels%theta(1:n_parcels), start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, qv_id, parcels%qv(1:n_parcels), start, cnt)
            call write_netcdf_dataset(ncid, ql_id, parcels%ql(1:n_parcels), start, cnt)
#endif
            call close_netcdf_file(ncid)

        end subroutine write_parcels

end program test_mpi_parcel_read_rejection
