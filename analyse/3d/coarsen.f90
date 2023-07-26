program coarsen_prog
    use constants, only : pi, twopi, f14, f12, zero, one, four
    use netcdf_reader
    use netcdf_writer
    implicit none

    ! Ordering in physical space: z, y, x

    ! Grid dimensions:
    integer :: fnx, fny, fnz
    integer :: cnx, cny, cnz

    ! Width and height of the domain:
    double precision :: extent(3), lower(3)
    integer          :: ncells(3)

    character(len=512) :: filename
    integer            :: shrink

!     call mpi_env_initialise

    call parse_command_line

    call get_domain

    if (mod(fnx, shrink) /= 0) then
        print '(i3, a18, i2, a1)', fnx, "not a multiple of", shrink, "."
        stop
    endif

    if (mod(fny, shrink) /= 0) then
        print '(i3, a18, i2, a1)', fny, "not a multiple of", shrink, "."
        stop
    endif

    if (mod(fnz, shrink) /= 0) then
        print '(i3, a18, i2, a1)', fnz, "not a multiple of", shrink, "."
        stop
    endif

    cnx = fnx / shrink
    cny = fny / shrink
    cnz = fnz / shrink

    print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', fnx, ' ny = ', fny, ' nz = ', fnz

    ! read data into array pp:
    call process_steps

    contains

        subroutine get_domain
            integer          :: ncid
            ! read domain dimensions
            call open_netcdf_file(trim(filename), NF90_NOWRITE, ncid)
            call get_netcdf_box(ncid, lower, extent, ncells)
            fnx = ncells(1)
            fny = ncells(2)
            fnz = ncells(3)
            call close_netcdf_file(ncid)
        end subroutine get_domain

        subroutine process_steps
            integer                       :: ncid, cnt(4), start(4), n_steps, ncerr
            integer                       :: mcid
            integer                       :: nDimensions, nVariables, step
            integer                       :: v, dimids(4), varids(32)
            character(len=64)             :: name
            character(len=64)             :: fields(32)
            double precision              :: fdata(fnz+1, fny, fnx)
            double precision              :: cdata(cnz+1, cny, cnx)
            character(512)                :: fname
            integer                       :: coord_ids(3)  ! = (x, y, z)
            integer                       :: t_axis_id
            double precision              :: dx(3)
            double precision, allocatable :: t(:)

            call open_netcdf_file(trim(filename), NF90_NOWRITE, ncid, l_serial=.true.)

            call get_num_steps(ncid, n_steps)

            allocate(t(n_steps))

            ncerr = nf90_inquire(ncid, nDimensions, nVariables)

            fname = 'coarsened_' // trim(filename)

            ncells(1) = cnx
            ncells(2) = cny
            ncells(3) = cnz

            call create_netcdf_file(ncfname=trim(fname),    &
                                    overwrite=.true.,       &
                                    ncid=mcid,              &
                                    l_serial=.true.)

            call write_netcdf_info(ncid=mcid,                    &
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)


            call write_netcdf_box(mcid, lower, extent, (/cnx, cny, cnz/))

!             call define_netcdf_dimension(mcid, "x", cnx,   dimids(1))
!             call define_netcdf_dimension(mcid, "y", cny,   dimids(2))
!             call define_netcdf_dimension(mcid, "z", cnz+1, dimids(3))

            call define_netcdf_spatial_dimensions_3d(ncid=mcid,                &
                                                     ngps=(/cnx, cny, cnz+1/), &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)

!             print *, "HI"
            call define_netcdf_temporal_dimension(mcid, dimids(4), t_axis_id)


            do v = 1, nVariables
                ncerr = nf90_inquire_variable(ncid, v, name)

                fields(v) = name
                if (trim(fields(v)) == 'x') then
                        cycle
                    endif

                    if (trim(fields(v)) == 'y') then
                        cycle
                    endif

                    if (trim(fields(v)) == 'z') then
                        cycle
                    endif

                    if (trim(fields(v)) == 't') then
                        cycle
                    endif

                print *, "Found " // trim(name)
                call define_netcdf_dataset(mcid,            &
                                           trim(name),      &
                                           '', '', '',      &
                                           NF90_DOUBLE,     &
                                           dimids,          &
                                           varids(v))
            enddo

            call close_definition(mcid)

            call open_netcdf_file(trim(fname), NF90_WRITE, mcid, l_serial=.true.)

            dx = extent / dble(ncells)
            call write_netcdf_axis_3d(mcid, dimids(1:3), lower, dx, &
                                      (/cnx, cny, cnz+1/))

            call read_netcdf_dataset(ncid, 't', t, start=(/1/), cnt=(/n_steps/))

            do step = 1, n_steps
                start = (/ 1,  1,  1, step /)

                do v = 1, nVariables
                    cnt = (/ fnx, fny, fnz+1, 1/)

                    if (trim(fields(v)) == 'x') then
                        cycle
                    endif

                    if (trim(fields(v)) == 'y') then
                        cycle
                    endif

                    if (trim(fields(v)) == 'z') then
                        cycle
                    endif

                    if (trim(fields(v)) == 't') then
                        call write_netcdf_scalar(mcid, t_axis_id, t(step), step)
                        cycle
                    endif

                    call read_netcdf_dataset(ncid, trim(fields(v)), fdata, &
                                             start=start, cnt=cnt)

                    print *, step, fields(v)

                    call coarsen(fdata, cdata)

                    cnt = (/ cnx, cny, cnz+1, 1/)

                    call write_netcdf_dataset(mcid, v, cdata, start=start, cnt=cnt, l_serial=.true.)

                enddo
            enddo

            call close_netcdf_file(mcid, l_serial=.true.)
            call close_netcdf_file(ncid, l_serial=.true.)

            deallocate(t)

        end subroutine process_steps

        subroutine coarsen(fdata, cdata)
            double precision  :: fdata(fnz+1, fny, fnx)
            double precision  :: cdata(cnz+1, cny, cnx)
            integer           :: ix, iy, iz, jx, jy, jz
            double precision  :: fac

            fac = one / dble(shrink ** 3)

            do ix = 1, cnx - 1
                jx = shrink * (ix-1) + 1
                do iy = 1, cny - 1
                    jy = shrink * (iy - 1) + 1
                    do iz = 2, cnz
                        jz = shrink * (iz - 1) + 1
                        cdata(iz, iy, ix) = sum(fdata(jz:jz+shrink, jy:jy+shrink, jx:jx+shrink)) * fac
                    enddo
                enddo
            enddo

!             cdata(1,     :, :) = fdata(1,     :, :)
!             cdata(cnz+1, :, :) = fdata(fnz+1, :, :)

        end subroutine coarsen

        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer            :: i, stat
            character(len=512) :: arg
            logical            :: exists

            shrink = -1
            filename = ''
            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--filename') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    filename = trim(arg)
                else if (arg == '--shrink') then
                    i = i + 1
                    call get_command_argument(i, arg)

                    ! 1 October 2021
                    ! https://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90
                    read(arg, *, iostat=stat)  shrink
                    if (stat .ne. 0) then
                        print *, 'Error conversion failed.'
                        stop
                    endif
                else if (arg == '--help') then
                    print *, 'This program computes the power spectrum and writes it to file.'
                    print *, 'An EPIC field output must be provided to analyse.'
                    print *, 'Run code with "coarsen --filename [field file]" --shrink [coarsen shrink]'
                    stop
                endif
                i = i+1
            enddo

            if ((filename == '') .or. (shrink == -1) ) then
                print *, 'No file or step provided. Run code with "coarsen --help"'
                stop
            endif

            ! check if correct file is passed
            stat = index(trim(filename), '_fields.nc', back=.true.)
            if (stat == 0) then
                print *, "Error: No EPIC field output file provided."
                stop
            endif

            ! check if file exsits
            inquire(file=trim(filename), exist=exists)
            if (.not. exists) then
                print *, "Error: File '" // trim(filename) // "' does not exist."
                stop
            endif
        end subroutine parse_command_line

end program coarsen_prog
