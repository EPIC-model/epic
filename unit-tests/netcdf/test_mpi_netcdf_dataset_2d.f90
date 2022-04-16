! =============================================================================
!                       Test netCDF 2D dataset
!
!               This unit test checks to write 2D datasets to netCDF.
! =============================================================================
program test_netcdf_dataset_2d
    use unit_test
    use netcdf_writer
    use mpi_communicator
    implicit none

    integer, parameter :: nx = 2, ny = 5
    integer            :: ix, iy, ncid, dimids(2)
    integer            :: var_id1 = -1, var_id2 = -1, var_id3 = -1
    double precision, allocatable :: dset(:, :)
    logical            :: passed = .true.
    integer            :: xstart, xend, ystart, yend, xlen, ylen

    call mpi_comm_initialise

    xstart = mpi_rank * nx + 1
    xend   = (mpi_rank + 1) * nx
    ystart = 1
    yend   = ny

    xlen = xend - xstart + 1
    ylen = yend - ystart + 1

    allocate(dset(ystart:yend, xstart:xend))

    do ix = xstart, xend
        do iy = ystart, yend
            dset(iy, ix) = mpi_rank
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

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id1)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'y_velocity', '', '',  'm/s', NF90_DOUBLE, dimids, var_id2)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'nparg', '', '',  '1', NF90_INT, dimids, var_id3)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',    &
                          access_flag=NF90_WRITE, &
                          ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id1, dset, start = (/xstart, ystart/), cnt = (/xlen, ylen/))

    passed = (passed .and. (ncerr == 0))

    dset = 1.5d0 + dset

    call write_netcdf_dataset(ncid, var_id2, dset, start = (/xstart, ystart/), cnt = (/xlen, ylen/))

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id3, dset, start = (/xstart, ystart/), cnt = (/xlen, ylen/))

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LOR, mpi_master, comm, mpi_err)
    endif

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI netCDF write 2D dataset', passed)
    endif

    deallocate(dset)

    call mpi_comm_finalise

end program test_netcdf_dataset_2d
