! =============================================================================
!          Test parallel parcel reading using the rejection method.
!
!   If a parcel file does not contain the start_index dataset,
!   all MPI ranks read all parcels, but reject parcels that do not
!   belong to the sub-domain owned by the rank.
! =============================================================================
program test_mpi_parcel_read_rejection
    use unit_test
    use options, only : parcel, write_netcdf_options
    use constants, only : one, zero, f12
    use parcels_mod, only : parcels
    use parcel_netcdf
    use parcel_init, only : parcel_default, init_timer
    use mpi_environment
    use mpi_layout
    use netcdf_utils
    use netcdf_writer
    use physics, only : write_physical_quantities
    use config, only : package_version, cf_version
    use iomanip, only : zfill
    use parameters, only : lower, update_parameters, extent, nx, ny, nz
    use mpi_timer
    implicit none

    logical            :: passed = .true.
    double precision   :: res
    integer            :: n_parcels_before
    double precision   :: x_sum, y_sum, z_sum

    integer            :: ncid
    integer            :: npar_dim_id, vol_id, buo_id,      &
                          x_pos_id, y_pos_id, z_pos_id,     &
                          x_vor_id, y_vor_id, z_vor_id,     &
                          b11_id, b12_id, b13_id,           &
                          b22_id, b23_id,                   &
                          t_axis_id, t_dim_id, mpi_dim_id
#ifndef ENABLE_DRY_MODE
    integer :: hum_id
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

    parcel%n_per_cell = 8

    call update_parameters

    call register_timer('parcel I/O', parcel_io_timer)
    call register_timer('parcel init', init_timer)

    !
    ! write parcels first
    !

    call parcel_default

    n_parcels_before = parcels%local_num

    ! we cannot check individual positions since the order of the parcels can be arbitrary
    x_sum = sum(parcels%position(1, 1:parcels%local_num))
    y_sum = sum(parcels%position(2, 1:parcels%local_num))
    z_sum = sum(parcels%position(3, 1:parcels%local_num))

    parcels%B(:, 1:parcels%local_num) = world%rank + 1
    parcels%volume(1:parcels%local_num) = world%rank + 1
    parcels%vorticity(:, 1:parcels%local_num) = world%rank + 1
    parcels%buoyancy(1:parcels%local_num) = world%rank + 1

#ifndef ENABLE_DRY_MODE
    parcels%humidity(1:parcels%local_num) = world%rank + 1
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
    parcels%buoyancy = 0
#ifndef ENABLE_DRY_MODE
    parcels%humidity = 0
#endif

    call read_netcdf_parcels('nctest_0000000001_parcels.nc')

    passed = (passed .and. (parcels%local_num == n_parcels_before))

    if (passed) then
        res = dble(world%rank + 1)

        passed = (passed .and. (abs(sum(parcels%position(1, 1:parcels%local_num)) - x_sum) == zero))
        passed = (passed .and. (abs(sum(parcels%position(2, 1:parcels%local_num)) - y_sum) == zero))
        passed = (passed .and. (abs(sum(parcels%position(3, 1:parcels%local_num)) - z_sum) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(1, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(2, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(3, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(4, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%B(5, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%volume(1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(1, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(2, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%vorticity(3, 1:parcels%local_num) - res)) == zero))
        passed = (passed .and. (maxval(abs(parcels%buoyancy(1:parcels%local_num) - res)) == zero))
#ifndef ENABLE_DRY_MODE
        passed = (passed .and. (maxval(abs(parcels%humidity(1:parcels%local_num) - res)) == zero))
#endif
    endif

    call parcels%deallocate

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

            ! all cores must know the correct number of total parcels
            parcels%total_num = parcels%local_num
            call MPI_Allreduce(MPI_IN_PLACE,        &
                               parcels%total_num,   &
                               1,                   &
                               MPI_INTEGER,         &
                               MPI_SUM,             &
                               world%comm,          &
                               world%err)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                         &
                                         name='n_parcels',                  &
                                         dimsize=int(parcels%total_num),    &
                                         dimid=npar_dim_id)

            call define_netcdf_dimension(ncid=ncid,                   &
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
                                       name='buoyancy',                         &
                                       long_name='parcel buoyancy',             &
                                       std_name='',                             &
                                       unit='m/s^2',                            &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=buo_id)

#ifndef ENABLE_DRY_MODE
            call define_netcdf_dataset(ncid=ncid,                               &
                                       name='humidity',                         &
                                       long_name='parcel humidity',             &
                                       std_name='',                             &
                                       unit='1',                                &
                                       dtype=NF90_DOUBLE,                       &
                                       dimids=dimids,                           &
                                       varid=hum_id)
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
            sendbuf(world%rank+1:world%size) = parcels%local_num
            sendbuf(world%rank+1) = 0

            call MPI_Reduce_scatter(sendbuf, start_index, recvcounts, MPI_INT, MPI_SUM, world%comm, world%err)

            ! we need to increase the start_index by 1
            ! since the starting index in Fortran is 1 and not 0.
            start_index = start_index + 1

            start = (/ start_index,       1 /)
            cnt   = (/ parcels%local_num, 1 /)

            call write_netcdf_dataset(ncid, x_pos_id, parcels%position(1, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, y_pos_id, parcels%position(2, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, z_pos_id, parcels%position(3, 1:parcels%local_num), start, cnt)

            call write_netcdf_dataset(ncid, b11_id, parcels%B(1, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, b12_id, parcels%B(2, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, b13_id, parcels%B(3, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, b22_id, parcels%B(4, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, b23_id, parcels%B(5, 1:parcels%local_num), start, cnt)

            call write_netcdf_dataset(ncid, vol_id, parcels%volume(1:parcels%local_num), start, cnt)

            call write_netcdf_dataset(ncid, x_vor_id, parcels%vorticity(1, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, y_vor_id, parcels%vorticity(2, 1:parcels%local_num), start, cnt)
            call write_netcdf_dataset(ncid, z_vor_id, parcels%vorticity(3, 1:parcels%local_num), start, cnt)

            call write_netcdf_dataset(ncid, buo_id, parcels%buoyancy(1:parcels%local_num), start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_netcdf_dataset(ncid, hum_id, parcels%humidity(1:parcels%local_num), start, cnt)
#endif
            call close_netcdf_file(ncid)

        end subroutine write_parcels

end program test_mpi_parcel_read_rejection
