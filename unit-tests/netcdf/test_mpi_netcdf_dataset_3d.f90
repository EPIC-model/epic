! =============================================================================
!                       Test MPI netCDF 3D dataset
!
!        This unit test checks to write 3D datasets in parallel to netCDF.
! =============================================================================
program test_mpi_netcdf_dataset_3d
    use unit_test
    use netcdf_writer
    use mpi_communicator, only : comm_world, mpi_err, mpi_comm_initialise, mpi_comm_finalise, comm_cart
    implicit none

    integer, parameter   :: nx = 10, ny = 8, nz = 5
    integer              :: ix, iy, iz, ncid, dimids(3)
    integer              :: var_id1 = -1, var_id2 = -1, var_id3
    integer, allocatable :: dset(:, :, :)
    logical              :: passed = .true.
    integer              :: start(3), cnt(3)
    integer              :: dims(2), coords(2)
    logical              :: periods(2)
    integer              :: rank ! we do not reorder the rank numbers, so this is unused!
    integer              :: lo(3), hi(3)

    call mpi_comm_initialise

    ! create slabs, z-direction keeps 1 processor
    dims = (/0, 0/)
    call MPI_Dims_create(mpi_size, 2, dims, mpi_err)

    periods = (/.true., .true./)

    call MPI_Cart_create(comm_world, 2, dims, periods, .false., comm_cart, mpi_err)

    call MPI_Comm_rank(comm_cart, rank, mpi_err)

    call MPI_Cart_coords(comm_cart, rank, 2, coords)

    call set_local_bounds(nx, coords(1), dims(1), lo(1), hi(1))
    call set_local_bounds(ny, coords(2), dims(2), lo(2), hi(2))
    lo(3) = 0
    hi(3) = nz

    allocate(dset(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)))

    do ix = lo(1), hi(1)
        do iy = lo(2), hi(2)
            do iz = lo(3), hi(3)
                dset(iz, iy, ix) = mpi_rank
            enddo
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc',                    &
                            overwrite=.true.,                       &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "x", nx, dimids(1))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "y", ny, dimids(2))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "z", nz+1, dimids(3))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id1)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'y_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id2)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'nparg', '', '', '-', NF90_INT, dimids, var_id3)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call open_netcdf_file(ncfname='nctest.nc',          &
                          access_flag=NF90_WRITE,       &
                          ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    start = lo + 1      ! need to add 1 since start must begin with index 1
    cnt = hi - lo + 1

    call write_netcdf_dataset(ncid, var_id1, dset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    dset = 100 + dset

    call write_netcdf_dataset(ncid, var_id2, dset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call write_netcdf_dataset(ncid, var_id3, dset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (mpi_rank == mpi_master) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, mpi_master, comm_world, mpi_err)
    endif

    if (mpi_rank == mpi_master) then
        call print_result_logical('Test MPI netCDF write 3D datasets', passed)
    endif

    deallocate(dset)

    call mpi_comm_finalise

    contains

        subroutine set_local_bounds(nglobal, coords, dims, first, last)
            integer, intent(in)  :: nglobal, coords, dims
            integer, intent(out) :: first, last
            integer              :: nlocal, remaining

            nlocal = nglobal / dims
            remaining = nglobal - dims * nlocal
            first = nlocal * coords
            if (coords < remaining) then
                nlocal = nlocal + 1
                first = first + coords
                last = last + coords
            else
                first = first + remaining
            endif
            last = first + nlocal - 1
        end subroutine set_local_bounds

end program test_mpi_netcdf_dataset_3d
