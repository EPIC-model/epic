! =============================================================================
!                       Test MPI netCDF reading 3D datasets
!
!             This unit test checks to read 3D datasets in parallel.
! =============================================================================
program test_mpi_netcdf_read_dataset_3d
    use unit_test
    use netcdf_writer
    use netcdf_reader
    use mpi_environment, only : world, mpi_env_initialise, mpi_env_finalise
    use mpi_layout, only : cart
    implicit none

    integer, parameter            :: nx = 5, ny = 10, nz = 15, nt = 3
    integer                       :: ix, iy, iz, ncid, dimids(4), t
    integer                       :: var_id = -1, cnt(4), start(4)
    double precision, allocatable :: wdset(:, :, :), rdset(:, :, :)
    double precision              :: error
    logical                       :: passed = .true.
    integer                       :: dims(2), coords(2)
    logical                       :: periods(2)
    integer                       :: lo(3), hi(3)

    call mpi_env_initialise

    ! create slabs, z-direction keeps 1 processor
    dims = (/0, 0/)
    call MPI_Dims_create(world%size, 2, dims, world%err)

    periods = (/.true., .true./)

    call MPI_Cart_create(world%comm, 2, dims, periods, .false., cart%comm, world%err)

    call MPI_Comm_rank(cart%comm, cart%rank, world%err)

    call MPI_Cart_coords(cart%comm, cart%rank, 2, coords)

    call set_local_bounds(nx, coords(1), dims(1), lo(1), hi(1))
    call set_local_bounds(ny, coords(2), dims(2), lo(2), hi(2))
    lo(3) = 0
    hi(3) = nz

    allocate(wdset(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)))
    allocate(rdset(lo(3):hi(3), lo(2):hi(2), lo(1):hi(1)))

    do ix = lo(1), hi(1)
        do iy = lo(2), hi(2)
            do iz = lo(3), hi(3)
                wdset(iz, iy, ix) = iy + ix * ny + 3 * iz
            enddo
        enddo
    enddo

    call create_netcdf_file(ncfname='nctest.nc', &
                            overwrite=.true.,    &
                            ncid=ncid)

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "x", nx, dimids(1))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "y", ny, dimids(2))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "z", nz+1, dimids(3))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dimension(ncid, "t", NF90_UNLIMITED, dimids(4))

    passed = (passed .and. (ncerr == 0))

    call define_netcdf_dataset(ncid, 'x_velocity', '', '', 'm/s', NF90_DOUBLE, dimids, var_id)

    passed = (passed .and. (ncerr == 0))

    call close_definition(ncid)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    ! write data

    do t = 1, nt

        start(1:3) = lo + 1      ! need to add 1 since start must begin with index 1
        cnt(1:3) = hi - lo + 1

        cnt(4) = 1
        start(4) = t

        call open_netcdf_file(ncfname='nctest.nc',    &
                            access_flag=NF90_WRITE, &
                            ncid=ncid)

        passed = (passed .and. (ncerr == 0))

        call write_netcdf_dataset(ncid, var_id, wdset, start, cnt)

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

    start(1:3) = lo + 1      ! need to add 1 since start must begin with index 1
    cnt(1:3) = hi - lo + 1
    cnt(4) = 1
    start(4) = nt

    call read_netcdf_dataset(ncid, 'x_velocity', rdset, start, cnt)

    passed = (passed .and. (ncerr == 0))

    call close_netcdf_file(ncid)

    passed = (passed .and. (ncerr == 0))

    error = sum(abs(wdset - rdset))

    passed = (passed .and. (ncerr == 0) .and. (error == 0.0d0))

    call delete_netcdf_file(ncfname='nctest.nc')

    passed = (passed .and. (ncerr == 0))

    if (world%rank == world%root) then
        call MPI_Reduce(MPI_IN_PLACE, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    else
        call MPI_Reduce(passed, passed, 1, MPI_LOGICAL, MPI_LAND, world%root, world%comm, world%err)
    endif

    if (world%rank == world%root) then
        call print_result_logical('Test MPI netCDF read 3D dataset', passed)
    endif

    deallocate(wdset)
    deallocate(rdset)

    call mpi_env_finalise

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

end program test_mpi_netcdf_read_dataset_3d
