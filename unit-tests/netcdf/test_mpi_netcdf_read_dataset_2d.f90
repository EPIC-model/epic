! =============================================================================
!                       Test MPI netCDF reading 2D datasets
!
!               This unit test checks to read 2D datasets in parallel.
! =============================================================================
program test_mpi_netcdf_read_dataset_2d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    use mpi_communicator
    implicit none

    integer, parameter            :: nx = 5, ny = 10, nt = 3
    integer                       :: ix, iy, ncid, dimids(3), t
    integer                       :: var_id = -1, cnt(3), start(3)
    double precision, allocatable :: wdset(:, :), rdset(:, :), error
    logical                       :: passed = .true.
    integer                       :: xstart, xend, ystart, yend, xlen, ylen

    call mpi_comm_initialise

    xstart = mpi_rank * nx + 1
    xend   = (mpi_rank + 1) * nx
    ystart = 1
    yend   = ny

    xlen = xend - xstart + 1
    ylen = yend - ystart + 1

    allocate(wdset(ystart:yend, xstart:xend))
    allocate(rdset(ystart:yend, xstart:xend))

    do ix = xstart, xend
        do iy = ystart, yend
            wdset(iy, ix) = dble(mpi_rank + 1)
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "x", mpi_size * nx, dimids(1))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "y", ny, dimids(2))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, dimids(3))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))


    ! write data

    do t = 1, nt

        cnt   = (/ xlen,   ylen,    1 /)
        start = (/ xstart, ystart,  t /)

        call open_netcdf_file(ncfname='nctest.nc',    &
                            access_flag=NF90_WRITE, &
                            ncid=ncid)

        passed = (passed .and. (ncerr == 0))

        call write_netcdf_dataset(ncid, var_id, wdset, start=start, cnt=cnt)

        passed = (passed .and. (ncerr == 0))

        call close_netcdf_file(ncid)

        passed = (passed .and. (ncerr == 0))

        wdset = 1.0d0 + wdset
    enddo

    ! we need to subtract the last increase again
    ! since it was not written
    wdset = wdset - 1.0d0

    ! read last data
    call open_netcdf_file(ncfname='nctest.nc',        &
                            access_flag=NF90_NOWRITE, &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    cnt   = (/ xlen,   ylen,    1  /)
    start = (/ xstart, ystart,  nt /)

    call read_netcdf_dataset(ncid, 'x_velocity', rdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    error = sum(abs(wdset - rdset))

    passed = (passed .and. (ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    endif

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI netCDF read 2D dataset', passed)
    endif

    deallocate(wdset)
    deallocate(rdset)

    call mpi_comm_finalise

end program test_mpi_netcdf_read_dataset_2d
