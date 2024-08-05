module parcel_netcdf
    use options, only : output, verbose
    use constants, only : one
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader
    use parcel_container, only : parcels            &
                               , n_parcels          &
                               , n_total_parcels    &
                               , parcel_delete
    use parameters, only : nx, ny, nz, extent, lower, max_num_parcels, write_zeta_boundary_flag
    use config, only : package_version, cf_version
    use mpi_timer, only : start_timer, stop_timer
    use iomanip, only : zfill
    use options, only : write_netcdf_options
    use physics, only : write_physical_quantities
    use mpi_environment
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_layout, only : box
    use mpi_ops, only : MPI_SUM_64BIT
    use datatypes, only : int64
    use mpi_utils, only : mpi_exit_on_error, mpi_print, mpi_check_for_error
    use fields, only : is_contained
    implicit none

    private

    integer :: parcel_io_timer

    integer :: n_writes = 1
    character(len=512) :: ncbasename

    character(len=512) :: ncfname
    integer            :: ncid
    integer            :: npar_dim_id   &
                        , t_axis_id     &
                        , t_dim_id      &
                        , mpi_dim_id
    double precision   :: restart_time

    integer, parameter :: NC_START = 1  &
                        , NC_XLO   = 2  &
                        , NC_XHI   = 3  &
                        , NC_YLO   = 4  &
                        , NC_YHI   = 5  &
                        , NC_VOL   = 6  &
                        , NC_X_POS = 7  &
                        , NC_Y_POS = 8  &
                        , NC_Z_POS = 9  &
                        , NC_BUOY  = 10 &
                        , NC_X_VOR = 11 &
                        , NC_Y_VOR = 12 &
                        , NC_Z_VOR = 13 &
                        , NC_B11   = 14 &
                        , NC_B12   = 15 &
                        , NC_B13   = 16 &
                        , NC_B22   = 17 &
                        , NC_B23   = 18 &
                        , NC_B33   = 19

    logical :: l_first_write = .true.
    logical :: l_unable = .false.

#ifndef ENABLE_DRY_MODE
    integer, parameter :: NC_HUM   = 20

#ifdef ENABLE_LABELS
    integer, parameter :: NC_LABEL    = 21 &
                       ,  NC_DILUTION = 22

    type(netcdf_info) :: nc_dset(NC_DILUTION)
#else
    type(netcdf_info) :: nc_dset(NC_HUM)
#endif

#else

#ifdef ENABLE_LABELS
    integer, parameter :: NC_LABEL    = 20 &
                       ,  NC_DILUTION = 21

    type(netcdf_info) :: nc_dset(NC_DILUTION)
#else
    type(netcdf_info) :: nc_dset(NC_B33)
#endif
#endif

    public :: create_netcdf_parcel_file &
            , write_netcdf_parcels      &
            , read_netcdf_parcels       &
            , parcel_io_timer

    contains

        ! Create the parcel file.
        ! @param[in] basename of the file
        ! @param[in] overwrite the file
        subroutine create_netcdf_parcel_file(basename, overwrite, l_restart)
            character(*), intent(in)  :: basename
            logical,      intent(in)  :: overwrite
            logical,      intent(in)  :: l_restart
            logical                   :: l_exist
            integer                   :: dimids(2)
            integer                   :: n, n_total

            call set_netcdf_parcel_output

            ncfname =  basename // '_' // zfill(n_writes) // '_parcels.nc'

            ncbasename = basename

            restart_time = -one

            if (l_restart) then
                ! find the last parcel file in order to set "n_writes" properly
                call exist_netcdf_file(ncfname, l_exist)
                do while (l_exist)
                    n_writes = n_writes + 1
                    ncfname =  basename // '_' // zfill(n_writes) // '_parcels.nc'
                    call exist_netcdf_file(ncfname, l_exist)
                    if (l_exist) then
                        call open_netcdf_file(ncfname, NF90_NOWRITE, ncid)
                        call get_time(ncid, restart_time)
                        call close_netcdf_file(ncid)
                    endif
                enddo
                return
            endif

            ! all cores must know the correct number of total parcels
            n_total_parcels = n_parcels
            call MPI_Allreduce(MPI_IN_PLACE, n_total_parcels, 1, MPI_INTEGER_64BIT, &
                               MPI_SUM_64BIT, world%comm, world%err)

            if ((world%rank == world%root) .and. (n_total_parcels > huge(n_total))) then
                print *, "WARNING: Unable to write parcels to the NetCDF file"
                print *, "         as the number of total parcel exceeds integer limit."
                l_unable = .true.
                return
            endif
            l_unable = .false.

            n_total = int(n_total_parcels)

            call create_netcdf_file(ncfname, overwrite, ncid)

            ! define global attributes
            call write_netcdf_info(ncid=ncid,                    &
                                   version_tag=package_version,  &
                                   file_type='parcels',          &
                                   cf_version=cf_version)

            call write_netcdf_box(ncid, lower, extent, (/nx, ny, nz/))

            call write_physical_quantities(ncid)

            call write_netcdf_options(ncid)

            ! define dimensions
            call define_netcdf_dimension(ncid=ncid,                 &
                                         name='n_parcels',          &
                                         dimsize=n_total,           &
                                         dimid=npar_dim_id)

            call define_netcdf_dimension(ncid=ncid,                  &
                                         name='world%size',          &
                                         dimsize=world%size,         &
                                         dimid=mpi_dim_id)

            call define_netcdf_temporal_dimension(ncid, t_dim_id, t_axis_id)

            dimids = (/npar_dim_id, t_dim_id/)

            ! define parcel attributes
            do n = NC_START, NC_YHI
                if (nc_dset(n)%l_enabled) then
                    call define_netcdf_dataset(ncid=ncid,                       &
                                               name=nc_dset(n)%name,            &
                                               long_name=nc_dset(n)%long_name,  &
                                               std_name=nc_dset(n)%std_name,    &
                                               unit=nc_dset(n)%unit,            &
                                               dtype=nc_dset(n)%dtype,          &
                                               dimids=(/mpi_dim_id/),           &
                                               varid=nc_dset(n)%varid)
                endif
            enddo

            do n = NC_VOL, size(nc_dset)
                if (nc_dset(n)%l_enabled) then
                    call define_netcdf_dataset(ncid=ncid,                       &
                                               name=nc_dset(n)%name,            &
                                               long_name=nc_dset(n)%long_name,  &
                                               std_name=nc_dset(n)%std_name,    &
                                               unit=nc_dset(n)%unit,            &
                                               dtype=nc_dset(n)%dtype,          &
                                               dimids=dimids,                   &
                                               varid=nc_dset(n)%varid)
                endif
            enddo

            call close_definition(ncid)

            call close_netcdf_file(ncid)

        end subroutine create_netcdf_parcel_file

        ! Write parcels of the current time step into the parcel file.
        ! @param[in] t is the time
        subroutine write_netcdf_parcels(t)
            double precision, intent(in) :: t
            integer                      :: cnt(2), start(2)
            integer                      :: recvcounts(world%size)
            integer                      :: sendbuf(world%size), start_index
#ifdef ENABLE_LABELS
            integer                      :: n
#endif

            call start_timer(parcel_io_timer)

            if (t <= restart_time) then
                call stop_timer(parcel_io_timer)
                return
            endif

            call create_netcdf_parcel_file(trim(ncbasename), .true., .false.)

            if (l_unable) then
                call stop_timer(parcel_io_timer)
                return
            endif

            call open_netcdf_file(ncfname, NF90_WRITE, ncid)

            ! we must write the boundary flag here
            call write_zeta_boundary_flag(ncid)

            ! write time
            call write_netcdf_scalar(ncid, t_axis_id, t, 1)

            ! after this operation all MPI ranks know their starting index
            recvcounts = 1
            start_index = 0
            sendbuf = 0
            sendbuf(world%rank+1:world%size) = n_parcels
            sendbuf(world%rank+1) = 0

            call MPI_Reduce_scatter(sendbuf(1:world%size),  &
                                    start_index,            &
                                    recvcounts,             &
                                    MPI_INTEGER,            &
                                    MPI_SUM,                &
                                    world%comm,             &
                                    world%err)

            call mpi_check_for_error(world, &
                "in MPI_Reduce_scatter of parcel_netcdf::write_netcdf_parcels.")

            ! we need to increase the start_index by 1
            ! since the starting index in Fortran is 1 and not 0.
            start_index = start_index + 1

            start = (/ start_index, 1 /)
            cnt   = (/ n_parcels,   1 /)

            if (nc_dset(NC_START)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_START)%varid, (/start_index/), &
                                          start=(/1+world%rank, 1/), cnt=(/1, 1/))
            endif

            if (nc_dset(NC_XLO)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_XLO)%varid, (/box%lo(1)/), &
                                          start=(/1+world%rank, 1/), cnt=(/1, 1/))
            endif

            if (nc_dset(NC_XHI)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_XHI)%varid, (/box%hi(1)/), &
                                          start=(/1+world%rank, 1/), cnt=(/1, 1/))
            endif

            if (nc_dset(NC_YLO)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_YLO)%varid, (/box%lo(2)/), &
                                          start=(/1+world%rank, 1/), cnt=(/1, 1/))
            endif

            if (nc_dset(NC_YHI)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(NC_YHI)%varid, (/box%hi(2)/), &
                                          start=(/1+world%rank, 1/), cnt=(/1, 1/))
            endif

            call write_parcel_attribute(NC_X_POS, parcels%position(1, :), start, cnt)
            call write_parcel_attribute(NC_Y_POS, parcels%position(2, :), start, cnt)
            call write_parcel_attribute(NC_Z_POS, parcels%position(3, :), start, cnt)

            call write_parcel_attribute(NC_B11, parcels%B(1, :), start, cnt)
            call write_parcel_attribute(NC_B12, parcels%B(2, :), start, cnt)
            call write_parcel_attribute(NC_B13, parcels%B(3, :), start, cnt)
            call write_parcel_attribute(NC_B22, parcels%B(4, :), start, cnt)
            call write_parcel_attribute(NC_B23, parcels%B(5, :), start, cnt)
            call write_parcel_attribute(NC_B33, parcels%B(6, :), start, cnt)

            call write_parcel_attribute(NC_VOL, parcels%volume, start, cnt)

            call write_parcel_attribute(NC_X_VOR, parcels%vorticity(1, :), start, cnt)
            call write_parcel_attribute(NC_Y_VOR, parcels%vorticity(2, :), start, cnt)
            call write_parcel_attribute(NC_Z_VOR, parcels%vorticity(3, :), start, cnt)

            call write_parcel_attribute(NC_BUOY, parcels%buoyancy, start, cnt)

#ifndef ENABLE_DRY_MODE
            call write_parcel_attribute(NC_HUM, parcels%humidity, start, cnt)
#endif
#ifdef ENABLE_LABELS
            call write_parcel_attribute_int(NC_LABEL, parcels%label, start, cnt)
            call write_parcel_attribute(NC_DILUTION, parcels%dilution, start, cnt)
            ! reset the labels to Fortran index which corresponds to current label
            ! reset the dilution to get this from time step to time step
            do n = 1, n_parcels
               parcels%label(n) = start_index + n - 1
               parcels%dilution(n) = 0
            end do
#endif
            ! increment counter
            n_writes = n_writes + 1

            call close_netcdf_file(ncid)

            call stop_timer(parcel_io_timer)

        end subroutine write_netcdf_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine write_parcel_attribute(id, pdata, start, cnt)
            integer,          intent(in) :: id
            double precision, intent(in) :: pdata(:)
            integer,          intent(in) :: cnt(2), start(2)

            if (nc_dset(id)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(id)%varid,  &
                                          pdata(1:n_parcels),       &
                                          start, cnt)
            endif
        end subroutine write_parcel_attribute

#ifdef ENABLE_LABELS
        subroutine write_parcel_attribute_int(id, pdata, start, cnt)
            integer,          intent(in) :: id
            integer(kind=8),  intent(in) :: pdata(:)
            integer,          intent(in) :: cnt(2), start(2)

            if (nc_dset(id)%l_enabled) then
                call write_netcdf_dataset(ncid, nc_dset(id)%varid,  &
                                          pdata(1:n_parcels),       &
                                          start, cnt)
            endif
        end subroutine write_parcel_attribute_int
#endif

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine read_netcdf_parcels(fname)
            character(*),     intent(in) :: fname
            integer                      :: start_index, num_indices, end_index
            integer                      :: n, n_total, pfirst, plast
            integer                      :: avail_size, n_remaining, n_read
            integer                      :: start(2), xlo, xhi, ylo, yhi
            logical                      :: l_same_world_size

            call start_timer(parcel_io_timer)

            call set_netcdf_parcel_info

            call open_netcdf_file(fname, NF90_NOWRITE, ncid)

            call get_num_parcels(ncid, n_total)

            n_total_parcels = n_total

            ! The number of MPI ranks disagree! We cannot use the 'start_index'
            ! to read in parcels
            call get_dimension_size(ncid, 'world%size', num_indices)
            l_same_world_size = (num_indices == world%size)

            if (l_same_world_size .and. has_dataset(ncid, 'start_index')) then

                if (world%rank < world%size - 1) then
                    ! we must add +1 since the start index is 1
                    call read_netcdf_dataset(ncid, 'start_index', start, (/world%rank + 1/), (/2/))
                    start_index = start(1)
                    ! we must subtract 1, otherwise rank reads the first parcel of rank+1
                    end_index = start(2) - 1
                else
                    ! the last MPI rank must only read the start index
                    call read_netcdf_dataset(ncid, 'start_index', start_index, num_indices)
                    end_index = n_total
                endif

                n_parcels = end_index - start_index + 1

                if (n_parcels > max_num_parcels) then
                    print *, "Number of parcels exceeds limit of", max_num_parcels, &
                             ". You may increase parcel%size_factor. Exiting."
                    call mpi_exit_on_error
                endif

                call read_chunk(start_index, end_index, 1)
            else if (has_dataset(ncid, 'xlo') .and. has_dataset(ncid, 'xhi') .and. &
                     has_dataset(ncid, 'ylo') .and. has_dataset(ncid, 'yhi') .and. &
                     has_dataset(ncid, 'start_index')) then
                !
                ! READ PARCEL WITH REJECTION METHOD BUT MAKING USE OF
                ! MPI BOX LAYOUT
                !
                call mpi_print("WARNING: MPI ranks disagree. Reading parcels with optimised rejection method!")

                n_parcels = 0
                pfirst = 1

                do n = 1, num_indices

                    ! read local box
                    call read_netcdf_dataset(ncid, 'xlo', xlo, start=n)
                    call read_netcdf_dataset(ncid, 'xhi', xhi, start=n)
                    call read_netcdf_dataset(ncid, 'ylo', ylo, start=n)
                    call read_netcdf_dataset(ncid, 'yhi', yhi, start=n)

                    ! check if box overlap (19 April 2024, https://stackoverflow.com/a/3269471):
                    if ((xlo <= box%hi(1)) .and. (box%lo(1) <= xhi) .and. &
                        (ylo <= box%hi(2)) .and. (box%lo(2) <= yhi)) then

                        ! get start and end index:
                        if (n < num_indices) then
                            call read_netcdf_dataset(ncid, 'start_index', start, (/n/), (/2/))
                            start_index = start(1)
                            ! we must subtract 1, otherwise rank reads the first parcel of the next domain
                            end_index = start(2) - 1
                        else
                            ! for the last index we can only read the start index
                            call read_netcdf_dataset(ncid, 'start_index', start_index, num_indices)
                            call get_num_parcels(ncid, n_total)
                            end_index = n_total
                        endif

                        call rejection_method(start_index, end_index, pfirst)

                        ! set pfirst to the end of the parcel container
                        pfirst = n_parcels + 1
                    endif
                enddo
            else
                !
                ! READ PARCELS WITH REJECTION METHOD
                ! (reject all parcels that are not part of
                !  the sub-domain owned by *this* MPI rank)
                !
                call mpi_print("WARNING: Unable to retrieve information for fast parcel reading.")
                call mpi_print("         All MPI ranks read all parcels!")

                start_index = 1
                end_index = min(max_num_parcels, n_total)
                pfirst = 1
                n_remaining = n_total
                n_parcels = 0

                do while (start_index <= end_index)

                    n_read = end_index - start_index + 1
                    n_remaining = n_remaining - n_read

                    call rejection_method(start_index, end_index, pfirst)

                    ! adjust the chunk size to fit the remaining memory
                    ! in the parcel container
                    avail_size = max(0, max_num_parcels - n_parcels)

                    ! update start index to fill container
                    pfirst = 1 + n_parcels
                    plast = min(pfirst + avail_size, n_total, max_num_parcels)

                    ! we must make sure that we have enough data in the
                    ! file as well as in the parcel container
                    n_read = min(plast - pfirst, n_remaining)

                    ! update start and end index for reading chunk
                    start_index = 1 + end_index
                    end_index = end_index + n_read
                enddo
            endif

            call close_netcdf_file(ncid)

            ! verify result
            n_total = n_parcels
            call MPI_Allreduce(MPI_IN_PLACE,    &
                               n_total,         &
                               1,               &
                               MPI_INTEGER,     &
                               MPI_SUM,         &
                               world%comm,      &
                               world%err)

            call mpi_check_for_error(world, &
                "in MPI_Allreduce of parcel_netcdf::read_netcdf_parcels.")

            if (n_total_parcels .ne. n_total) then
                call mpi_exit_on_error(&
                    "Local number of parcels does not sum up to total number!")
            endif

            call stop_timer(parcel_io_timer)

        end subroutine read_netcdf_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine rejection_method(start_index, end_index, pfirst)
            integer, intent(in)  :: start_index
            integer, intent(in)  :: end_index
            integer, intent(in)  :: pfirst
            integer              :: m, k, n_read
            integer, allocatable :: invalid(:)

            call read_chunk(start_index, end_index, pfirst)
            n_read = end_index - start_index + 1
            n_parcels = n_parcels + n_read

            ! if all MPI ranks read all parcels, each MPI rank must delete the parcels
            ! not belonging to its domain
            allocate(invalid(0:n_read))

            m = 1
            do k = pfirst, n_parcels
                if (is_contained(parcels%position(:, k))) then
                    cycle
                endif

                invalid(m) = k

                m = m + 1
            enddo

            ! remove last increment
            m = m - 1

            ! updates the variable n_parcels
            call parcel_delete(invalid(0:m), n_del=m)

            deallocate(invalid)

        end subroutine rejection_method

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! This subroutine assumes the NetCDF file to be open.
        subroutine read_chunk(first, last, pfirst)
            integer, intent(in) :: first, last, pfirst
            logical             :: l_valid = .false.
            integer             :: cnt(2), start(2)
            integer             :: num, plast, n

            num = last - first + 1
            plast = pfirst + num - 1

            if (plast > max_num_parcels) then
                print *, "Number of parcels exceeds limit of", max_num_parcels, &
                         ". You may increase parcel%size_factor. Exiting."
                call mpi_exit_on_error
            endif

            start = (/ first, 1 /)
            cnt   = (/ num,   1 /)

            if (has_dataset(ncid, 'B11')) then
                call read_netcdf_dataset(ncid, 'B11', parcels%B(1, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B11 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'B12')) then
                call read_netcdf_dataset(ncid, 'B12', parcels%B(2, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B12 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'B13')) then
                call read_netcdf_dataset(ncid, 'B13', parcels%B(3, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B13 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'B22')) then
                call read_netcdf_dataset(ncid, 'B22', parcels%B(4, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B22 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'B23')) then
                call read_netcdf_dataset(ncid, 'B23', parcels%B(5, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B23 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'B33')) then
                call read_netcdf_dataset(ncid, 'B33', parcels%B(6, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel shape component B33 must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'x_position')) then
                call read_netcdf_dataset(ncid, 'x_position', &
                                         parcels%position(1, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel x position must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'y_position')) then
                call read_netcdf_dataset(ncid, 'y_position', &
                                         parcels%position(2, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel y position must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'z_position')) then
                call read_netcdf_dataset(ncid, 'z_position', &
                                         parcels%position(3, pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel z position must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'volume')) then
                call read_netcdf_dataset(ncid, 'volume', &
                                         parcels%volume(pfirst:plast), start, cnt)
            else
                call mpi_exit_on_error(&
                    "The parcel volume must be present! Exiting.")
            endif

            if (has_dataset(ncid, 'x_vorticity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'x_vorticity', &
                                         parcels%vorticity(1, pfirst:plast), start, cnt)
            endif

            if (has_dataset(ncid, 'y_vorticity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'y_vorticity', &
                                         parcels%vorticity(2, pfirst:plast), start, cnt)
            endif

            if (has_dataset(ncid, 'z_vorticity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'z_vorticity', &
                                         parcels%vorticity(3, pfirst:plast), start, cnt)
            endif

            if (has_dataset(ncid, 'buoyancy')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'buoyancy', &
                                         parcels%buoyancy(pfirst:plast), start, cnt)
            endif

#ifndef ENABLE_DRY_MODE
            if (has_dataset(ncid, 'humidity')) then
                l_valid = .true.
                call read_netcdf_dataset(ncid, 'humidity', &
                                         parcels%humidity(pfirst:plast), start, cnt)
            endif
#endif
#ifdef ENABLE_LABELS
            ! reset the labels to Fortran index which corresponds to current label
            ! reset the dilution to get this from time step to time step
            do n = pfirst, plast
               parcels%label(n) = first + n - pfirst 
               parcels%dilution(n) = 0
            end do
#endif


            if (.not. l_valid) then
                call mpi_exit_on_error(&
                    "Either the parcel buoyancy or vorticity must be present! Exiting.")
            endif
        end subroutine read_chunk

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_parcel_output
            integer :: n
            logical :: l_enabled_restart = .true.

            call set_netcdf_parcel_info

            ! check custom tags
            if (any('all' == output%parcel_list(:))) then
                nc_dset(:)%l_enabled = .true.
            else if (any('default' == output%parcel_list(:))) then
                nc_dset(NC_X_POS)%l_enabled = .true.
                nc_dset(NC_Y_POS)%l_enabled = .true.
                nc_dset(NC_Z_POS)%l_enabled = .true.
                nc_dset(NC_X_VOR)%l_enabled = .true.
                nc_dset(NC_Y_VOR)%l_enabled = .true.
                nc_dset(NC_Z_VOR)%l_enabled = .true.
                nc_dset(NC_BUOY)%l_enabled  = .true.
#ifndef ENABLE_DRY_MODE
                nc_dset(NC_HUM)%l_enabled   = .true.
#endif
#ifdef ENABLE_LABELS
                nc_dset(NC_LABEL)%l_enabled   = .true.
                nc_dset(NC_DILUTION)%l_enabled   = .true.
#endif
                nc_dset(NC_VOL)%l_enabled   = .true.
                nc_dset(NC_START)%l_enabled = .true.
                nc_dset(NC_XLO)%l_enabled   = .true.
                nc_dset(NC_XHI)%l_enabled   = .true.
                nc_dset(NC_YLO)%l_enabled   = .true.
                nc_dset(NC_YHI)%l_enabled   = .true.
                nc_dset(NC_B11)%l_enabled   = .true.
                nc_dset(NC_B12)%l_enabled   = .true.
                nc_dset(NC_B13)%l_enabled   = .true.
                nc_dset(NC_B22)%l_enabled   = .true.
                nc_dset(NC_B23)%l_enabled   = .true.
                nc_dset(NC_B33)%l_enabled   = .true.
            else
                ! check individual fields
                do n = 1, size(nc_dset)
                    nc_dset(n)%l_enabled = any(nc_dset(n)%name == output%parcel_list(:))
                enddo
            endif

            if (count(nc_dset(:)%l_enabled) == 0) then
                if ((world%rank == world%root) .and. l_first_write) then
                    print *, "WARNING: No parcel attributes are actively selected. EPIC is going to write"
                    print *, "         the default parcel attributes. Stop the simulation now if this is"
                    print *, "         not your intention. Parcel attributes can be provided to the list"
                    print *, "         'output%parcel_list' in the configuration file."
                    print *, "         The following parcel attributes are available:"
                    do n = 1, size(nc_dset)
                        print *, "         " // nc_dset(n)%name // " : " // trim(nc_dset(n)%long_name)
                    enddo
                    print *, "         " // "all"     // repeat(" ", 29) // " : write all parcel attributes"
                    print *, "         " // "default" // repeat(" ", 25) // " : write default parcel attributes"
                    print *, ""
                endif
                nc_dset(NC_START)%l_enabled = .true.
                nc_dset(NC_XLO)%l_enabled   = .true.
                nc_dset(NC_XHI)%l_enabled   = .true.
                nc_dset(NC_YLO)%l_enabled   = .true.
                nc_dset(NC_YHI)%l_enabled   = .true.
                nc_dset(NC_X_POS)%l_enabled = .true.
                nc_dset(NC_Y_POS)%l_enabled = .true.
                nc_dset(NC_Z_POS)%l_enabled = .true.
                nc_dset(NC_X_VOR)%l_enabled = .true.
                nc_dset(NC_Y_VOR)%l_enabled = .true.
                nc_dset(NC_Z_VOR)%l_enabled = .true.
                nc_dset(NC_BUOY)%l_enabled  = .true.
#ifndef ENABLE_DRY_MODE
                nc_dset(NC_HUM)%l_enabled   = .true.
#endif
#ifdef ENABLE_LABELS
                nc_dset(NC_LABEL)%l_enabled   = .true.
                nc_dset(NC_DILUTION)%l_enabled   = .true.
#endif
                nc_dset(NC_VOL)%l_enabled   = .true.
                nc_dset(NC_B11)%l_enabled   = .true.
                nc_dset(NC_B12)%l_enabled   = .true.
                nc_dset(NC_B13)%l_enabled   = .true.
                nc_dset(NC_B22)%l_enabled   = .true.
                nc_dset(NC_B23)%l_enabled   = .true.
                nc_dset(NC_B33)%l_enabled   = .true.
            endif

#ifdef ENABLE_VERBOSE
            if (verbose .and. (world%rank == world%root) .and. l_first_write) then
                print *, "EPIC is going to write the following parcel attributes:"
                do n = 1, size(nc_dset)
                    if (nc_dset(n)%l_enabled) then
                        print *, repeat(" ", 4) // trim(nc_dset(n)%name)
                    endif
                enddo
                print *, ""
            endif
#endif


            if (l_first_write) then
                l_first_write = .false.

                do n = NC_X_POS, NC_Y_POS
                    l_enabled_restart = (l_enabled_restart .and. nc_dset(n)%l_enabled)
                enddo

                do n = NC_B11, NC_B33
                    l_enabled_restart = (l_enabled_restart .and. nc_dset(n)%l_enabled)
                enddo
                l_enabled_restart = (l_enabled_restart .and. nc_dset(NC_VOL)%l_enabled)

                l_enabled_restart = (l_enabled_restart .and.                &
                                      ((nc_dset(NC_X_VOR)%l_enabled .and.   &
                                        nc_dset(NC_Y_VOR)%l_enabled .and.   &
                                        nc_dset(NC_Z_VOR)%l_enabled) .or.   &
                                        nc_dset(NC_BUOY)%l_enabled))


                if ((.not. l_enabled_restart) .and. (world%rank == world%root)) then
                    print *, "WARNING: EPIC will not be able to restart from a parcel file."
                    print *, "         You must at least write the B-shape matrix, parcel position"
                    print *, "         parcel volume and parcel vorticity or buoyancy to enable a"
                    print *, "         restart. If you intend to restart from a parcel file later,"
                    print *, "         you must stop the simulation immediately. Furthermore, you can"
                    print *, "         write the MPI 'start_index' to speed up the restart."
                endif
            endif

        end subroutine set_netcdf_parcel_output

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine set_netcdf_parcel_info

            nc_dset(NC_X_POS) = netcdf_info(name='x_position',                   &
                                            long_name='x position component',    &
                                            std_name='',                         &
                                            unit='m',                            &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_Y_POS) = netcdf_info(name='y_position',                   &
                                            long_name='y position component',    &
                                            std_name='',                         &
                                            unit='m',                            &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_Z_POS) = netcdf_info(name='z_position',                   &
                                            long_name='z position component',    &
                                            std_name='',                         &
                                            unit='m',                            &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_START) = netcdf_info(name='start_index',                  &
                                            long_name='MPI rank start index',    &
                                            std_name='',                         &
                                            unit='1',                            &
                                            dtype=NF90_INT)

            nc_dset(NC_XLO) = netcdf_info(name='xlo',                           &
                                          long_name='lower box boundary in x',  &
                                          std_name='',                          &
                                          unit='1',                             &
                                          dtype=NF90_INT)

            nc_dset(NC_XHI) = netcdf_info(name='xhi',                           &
                                          long_name='upper box boundary in x',  &
                                          std_name='',                          &
                                          unit='1',                             &
                                          dtype=NF90_INT)

            nc_dset(NC_YLO) = netcdf_info(name='ylo',                           &
                                          long_name='lower box boundary in y',  &
                                          std_name='',                          &
                                          unit='1',                             &
                                          dtype=NF90_INT)

            nc_dset(NC_YHI) = netcdf_info(name='yhi',                           &
                                          long_name='upper box boundary in y',  &
                                          std_name='',                          &
                                          unit='1',                             &
                                          dtype=NF90_INT)

            nc_dset(NC_B11) = netcdf_info(name='B11',                              &
                                          long_name='B11 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_B12) = netcdf_info(name='B12',                              &
                                          long_name='B12 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_B13) = netcdf_info(name='B13',                              &
                                          long_name='B13 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_B22) = netcdf_info(name='B22',                              &
                                          long_name='B22 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_B23) = netcdf_info(name='B23',                              &
                                          long_name='B23 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_B33) = netcdf_info(name='B33',                              &
                                          long_name='B33 element of shape matrix', &
                                          std_name='',                             &
                                          unit='m^2',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_VOL) = netcdf_info(name='volume',                           &
                                          long_name='parcel volume',               &
                                          std_name='',                             &
                                          unit='m^3',                              &
                                          dtype=NF90_DOUBLE)

            nc_dset(NC_X_VOR) = netcdf_info(name='x_vorticity',                      &
                                            long_name='x vorticity component',       &
                                            std_name='',                             &
                                            unit='1/s',                              &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_Y_VOR) = netcdf_info(name='y_vorticity',                      &
                                            long_name='y vorticity component',       &
                                            std_name='',                             &
                                            unit='1/s',                              &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_Z_VOR) = netcdf_info(name='z_vorticity',                      &
                                            long_name='z vorticity component',       &
                                            std_name='',                             &
                                            unit='1/s',                              &
                                            dtype=NF90_DOUBLE)

            nc_dset(NC_BUOY) = netcdf_info(name='buoyancy',                         &
                                           long_name='parcel buoyancy',             &
                                           std_name='',                             &
                                           unit='m/s^2',                            &
                                           dtype=NF90_DOUBLE)

#ifndef ENABLE_DRY_MODE
            nc_dset(NC_HUM) = netcdf_info(name='humidity',                         &
                                          long_name='parcel humidity',             &
                                          std_name='',                             &
                                          unit='1',                                &
                                          dtype=NF90_DOUBLE)
#endif
#ifdef ENABLE_LABELS
            nc_dset(NC_LABEL) = netcdf_info(name='label',                          &
                                          long_name='parcel label',                &
                                          std_name='',                             &
                                          unit='1',                                &
                                          dtype=NF90_INT64)

            nc_dset(NC_DILUTION) = netcdf_info(name='dilution',                    &
                                          long_name='parcel log dilution',         &
                                          std_name='',                             &
                                          unit='1',                                &
                                          dtype=NF90_DOUBLE)
#endif
        end subroutine set_netcdf_parcel_info

end module parcel_netcdf
