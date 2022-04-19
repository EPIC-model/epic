! =============================================================================
! References:
! https://github.com/Unidata/netcdf-fortran/tree/master/examples (23.12.2021)
! https://docs.unidata.ucar.edu/netcdf-fortran/current/f90_The-NetCDF-Fortran-90-Interface-Guide.html
! =============================================================================
module netcdf_utils
    use netcdf
    use mpi_communicator
    implicit none

    ! netCDF error if non-zero
    integer :: ncerr = 0

    contains

        ! If the logical argument 'l_single = .true.', the file
        ! is created by the root MPI rank with MPI size > 1. The argument
        ! is .false. by default.
        subroutine create_netcdf_file(ncfname, overwrite, ncid, l_single)
            character(*),      intent(in)  :: ncfname
            logical,           intent(in)  :: overwrite
            integer,           intent(out) :: ncid
            logical                        :: l_exist = .true.
            logical, optional, intent(in)  :: l_single
            logical                        :: l_parallel = (mpi_size > 1)

            if (present(l_single)) then
                l_parallel = .not. l_single
            endif

            ! check whether file exists
            call exist_netcdf_file(ncfname, l_exist)

            if (l_exist .and. overwrite) then
                call delete_netcdf_file(ncfname)
            else if (l_exist) then
                print *, "File '" // trim(ncfname) // "' already exists. Exiting."
                stop
            endif

            if (l_parallel) then
                ! 16 April
                ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node25.htm
                ncerr = nf90_create(path = ncfname,                         &
                                    cmode = ior(NF90_NETCDF4, NF90_MPIIO),  &
                                    ncid = ncid,                            &
                                    comm = comm_world%MPI_VAL,              &
                                    info = MPI_INFO_NULL%MPI_VAL)
            else
                ! in single execution mpi_master = mpi_rank = 0
                if (mpi_master == mpi_rank) then
                    ncerr = nf90_create(path = ncfname,        &
                                        cmode = NF90_NETCDF4,  &
                                        ncid = ncid)
                endif
            endif

            call check_netcdf_error("Failed to create netcdf file'" // trim(ncfname) // "'.")

        end subroutine create_netcdf_file

        subroutine delete_netcdf_file(ncfname)
            character(*), intent(in) :: ncfname
            integer                  :: stat

            if (mpi_rank .ne. mpi_master) then
                return
            endif

            ! 15 June 2021
            ! https://stackoverflow.com/questions/18668832/how-delete-file-from-fortran-code
            open(unit=1234, iostat=stat, file=ncfname, status='old')
            if (stat == 0) then
                close(1234, status='delete')
            endif
        end subroutine delete_netcdf_file


        ! If the logical argument 'l_single = .true.', the file
        ! is opened by the root MPI rank with MPI size > 1. The argument
        ! is .false. by default.
        subroutine open_netcdf_file(ncfname, access_flag, ncid, l_single)
            character(*),      intent(in)  :: ncfname
            integer,           intent(in)  :: access_flag ! NF90_WRITE or NF90_NOWRITE
            integer,           intent(out) :: ncid
            logical, optional, intent(in)  :: l_single
            logical                        :: l_parallel = (mpi_size > 1)

            if (present(l_single)) then
                l_parallel = .not. l_single
            endif

            if (l_parallel) then
                ncerr = nf90_open(path = ncfname,               &
                                  mode = access_flag,           &
                                  ncid = ncid,                  &
                                  comm = comm_world%MPI_VAL,    &
                                  info = MPI_INFO_NULL%MPI_VAL)
            else
                ! in single execution mpi_master = mpi_rank = 0
                if (mpi_master == mpi_rank) then
                    ncerr = nf90_open(path = ncfname,     &
                                    mode = access_flag,   &
                                    ncid = ncid)
                endif
            endif

            call check_netcdf_error("Opening the netcdf file failed.")

        end subroutine open_netcdf_file

        ! Note: This subroutine should only be called with l_single = .true.
        ! if it was opened in single mode.
        subroutine close_netcdf_file(ncid, l_single)
            integer,           intent(in)  :: ncid
            logical, optional, intent(in)  :: l_single
            logical                        :: l_parallel = (mpi_size > 1)

            if (present(l_single)) then
                l_parallel = (.not. l_single)
            endif

            if (l_parallel) then
                ncerr = nf90_close(ncid)
            else
                if (mpi_master == mpi_rank) then
                    ncerr = nf90_close(ncid)
                endif
            endif

            call check_netcdf_error("Closing the netcdf file failed.")
        end subroutine close_netcdf_file

        subroutine exist_netcdf_file(ncfname, l_exist)
            character(*), intent(in)  :: ncfname
            logical,      intent(out) :: l_exist

            ! check whether file exists
            inquire(file=ncfname, exist=l_exist)
        end subroutine exist_netcdf_file

        subroutine check_netcdf_error(msg)
            character(*), intent(in) :: msg
#ifndef NDEBUG
            if (ncerr /= nf90_noerr .and. mpi_rank == mpi_master) then
                print *, msg
                print *, trim(nf90_strerror(ncerr))
                stop
            endif
#endif
        end subroutine check_netcdf_error
end module netcdf_utils
