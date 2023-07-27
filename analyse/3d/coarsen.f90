program coarsening
    use constants, only : f14, f12
    use netcdf_reader
    use netcdf_writer
    implicit none

    double precision, parameter :: f116 = 1.0d0 / 16.0d0
    double precision, parameter :: f38 = 3.0d0 / 8.0d0

    ! Grid dimensions:
    integer :: fnx, fny, fnz    ! fine grid
    integer :: cnx, cny, cnz    ! coarse grid

    ! Width and height of the domain:
    double precision :: extent(3), lower(3)
    integer          :: ncells(3)

    character(len=512) :: cname     ! file name of coarse grid
    character(len=512) :: fname     ! file name of fine grid
    integer            :: shrink

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

    select case(shrink)
        case (2)
            cname = 'crse_2_' // trim(fname)
        case (3)
            cname = 'crse_3_' // trim(fname)
        case default
            print *, "We can only coarse by a factor 2 or 3."
            stop
    end select

    print '(a23, i5, a6, i5, a6, i5)', 'Grid dimensions: nx = ', fnx, ' ny = ', fny, ' nz = ', fnz

    ! read data into array pp:
    call process_steps

    contains

        subroutine get_domain
            integer          :: ncid
            ! read domain dimensions
            call open_netcdf_file(trim(fname), NF90_NOWRITE, ncid)
            call get_netcdf_box(ncid, lower, extent, ncells)
            fnx = ncells(1)
            fny = ncells(2)
            fnz = ncells(3)
            call close_netcdf_file(ncid)
        end subroutine get_domain

        subroutine process_steps
            integer                       :: ncid, cnt(4), start(4), n_steps, ncerr
            integer                       :: mcid
            integer                       :: nDimensions, nVariables, step, nComp
            integer                       :: dimids(4), nc, nv
            character(len=64)             :: name
            integer, allocatable          :: varids(:)
            character(len=64), allocatable:: fields(:)
            double precision, allocatable :: fdata(:, :, :, :)
            double precision, allocatable :: cdata(:, :, :, :)
            integer                       :: coord_ids(3)  ! = (x, y, z)
            integer                       :: t_axis_id
            double precision              :: dx(3)
            double precision, allocatable :: t(:)

            call open_netcdf_file(trim(fname), NF90_NOWRITE, ncid, l_serial=.true.)

            call get_num_steps(ncid, n_steps)

            ncerr = nf90_inquire(ncid, nDimensions, nVariables)

            ! allocate for all variables (excluding spatial and temporal variables)
            nComp = nVariables - 4

            allocate(t(n_steps))
            allocate(fdata(0:fnz, 0:fny-1, 0:fnx-1, nComp))
            allocate(cdata(0:cnz, 0:cny-1, 0:cnx-1, nComp))
            allocate(fields(nComp))
            allocate(varids(nComp))

            ncells(1) = cnx
            ncells(2) = cny
            ncells(3) = cnz

            call create_netcdf_file(ncfname=trim(cname),    &
                                    overwrite=.true.,       &
                                    ncid=mcid,              &
                                    l_serial=.true.)

            call write_netcdf_info(ncid=mcid,                    &
                                   version_tag=package_version,  &
                                   file_type='fields',           &
                                   cf_version=cf_version)


            call write_netcdf_box(mcid, lower, extent, (/cnx, cny, cnz/))

            call define_netcdf_spatial_dimensions_3d(ncid=mcid,                &
                                                     ngps=(/cnx, cny, cnz+1/), &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)

            call define_netcdf_temporal_dimension(mcid, dimids(4), t_axis_id)


            !------------------------------------------------------------------
            ! Register all variables: We ignore the spatial and temporal
            ! dimensions.
            nc = 1
            do nv = 1, nVariables
                ncerr = nf90_inquire_variable(ncid, nv, name)

                select case(trim(name))
                    case ('x')
                        print *, "Skipping " // trim(name)
                    case ('y')
                        print *, "Skipping " // trim(name)
                    case ('z')
                        print *, "Skipping " // trim(name)
                    case ('t')
                        print *, "Skipping " // trim(name)
                    case default
                        print *, "Found " // trim(name)
                        fields(nc) = name
                        call define_netcdf_dataset(mcid,        &
                                                   trim(name),  &
                                                   '', '', '',  &
                                                   NF90_DOUBLE, &
                                                   dimids,      &
                                                   varids(nc))
                        nc = nc + 1
                end select
            enddo


            call close_definition(mcid)

            call close_netcdf_file(mcid, l_serial=.true.)

            call open_netcdf_file(trim(cname), NF90_WRITE, mcid, l_serial=.true.)

            dx = extent / dble(ncells)
            call write_netcdf_axis_3d(mcid,                 &
                                      dimids(1:3),          &
                                      lower,                &
                                      dx,                   &
                                      (/cnx, cny, cnz+1/))

            call read_netcdf_dataset(ncid, 't', t, start=(/1/), cnt=(/n_steps/))

            do step = 1, n_steps
                start = (/1, 1, 1, step/)

                !--------------------------------------------------------------
                ! read all fields:

                call write_netcdf_scalar(mcid, t_axis_id, t(step), step)

                cnt = (/fnx, fny, fnz+1, 1/)

                do nc = 1, nComp
                    call read_netcdf_dataset(ncid,                  &
                                             trim(fields(nc)),      &
                                             fdata(:, :, :, nc),    &
                                             start=start,           &
                                             cnt=cnt)
                enddo

                !--------------------------------------------------------------
                ! coarsen all fields:
                select case(shrink)
                    case (2)
                        call coarsen_two(fdata, cdata, nComp)
                    case (3)
                        call coarsen_three(fdata, cdata, nComp)
                    case default
                        print *, "We can only coarse by a factor 2 or 3."
                        stop
                end select

                !--------------------------------------------------------------
                ! write all coarsened fields:

                cnt = (/cnx, cny, cnz+1, 1/)

                do nc = 1, nComp
                    call write_netcdf_dataset(mcid,                 &
                                              nc,                   &
                                              cdata(:, :, :, nc),   &
                                              start=start,          &
                                              cnt=cnt,              &
                                              l_serial=.true.)
                enddo
            enddo

            call close_netcdf_file(mcid, l_serial=.true.)
            call close_netcdf_file(ncid, l_serial=.true.)

            deallocate(t)
            deallocate(fdata)
            deallocate(cdata)
            deallocate(fields)
            deallocate(varids)

        end subroutine process_steps

        subroutine coarsen_two(fdata, cdata, nc)
            integer, intent(in) :: nc
            double precision, intent(in)  :: fdata(0:fnz, 0:fny-1, 0:fnx-1, nc)
            double precision, intent(out) :: cdata(0:cnz, 0:cny-1, 0:cnx-1, nc)
            double precision              :: xdata(0:fnz, 0:fny-1, 0:cnx-1, nc)
            double precision              :: xydata(0:fnz, 0:cny-1, 0:cnx-1, nc)
            integer                       :: ix, jx, iy, iz, jy, jz, jxp1, jyp1, jxm1, jym1

            !------------------------------------------------------------------
            ! Coarsen in x using 1-2-1 stencil:
            do ix = 0, cnx-1
                jx = 2 * ix
                jxm1 = mod(jx - 1 + fnx, fnx)
                jxp1 = mod(jx + 1, fnx)

                do iy = 0, fny-1
                    do iz = 0, fnz
                        xdata(iz, iy, ix, :) = f14 * (fdata(iz, iy, jxm1, :) + &
                                                      fdata(iz, iy, jxp1, :))  &
                                             + f12 *  fdata(iz, iy, jx, :)
                    enddo
                enddo
            enddo

            !------------------------------------------------------------------
            ! Coarsen in y using 1-2-1 stencil:
            do ix = 0, cnx-1
                do iy = 0, cny-1
                    jy = 2 * iy
                    jym1 = mod(jy - 1 + fny, fny)
                    jyp1 = mod(jy + 1, fny)
                    do iz = 0, fnz
                        xydata(iz, iy, ix, :) = f14 * (xdata(iz, jym1, ix, :) + &
                                                       xdata(iz, jyp1, ix, :))  &
                                              + f12 *  xdata(iz, jy, ix, :)
                    enddo
                enddo
            enddo

            !------------------------------------------------------------------
            ! Coarsen in z using 1-2-1 stencil:
            ! Note: The boundaries are kept the same.
            do ix = 0, cnx-1
                do iy = 0, cny-1
                    do iz = 1, cnz - 1
                        jz = 2 * iz
                        cdata(iz, iy, ix, :) = f14 * (xydata(jz-1, iy, ix, :) + &
                                                      xydata(jz+1, iy, ix, :))  &
                                             + f12 *  xydata(jz, iy, ix, :)
                    enddo
                    cdata(0,   iy, ix, :) = xydata(0,   iy, ix, :)
                    cdata(cnz, iy, ix, :) = xydata(fnz, iy, ix, :)
                enddo
            enddo
        end subroutine coarsen_two

        subroutine coarsen_three(fdata, cdata, nc)
            integer, intent(in) :: nc
            double precision, intent(in)  :: fdata(0:fnz, 0:fny-1, 0:fnx-1, nc)
            double precision, intent(out) :: cdata(0:cnz, 0:cny-1, 0:cnx-1, nc)
            double precision              :: xdata(0:fnz, 0:fny-1, 0:cnx-1, nc)
            double precision              :: xydata(0:fnz, 0:cny-1, 0:cnx-1, nc)
            integer                       :: ix, jx, iy, iz, jy, jz, jxp1, jyp1, jxm1, jym1
            integer :: jxm2, jxp2
            integer :: jym2, jyp2

            !------------------------------------------------------------------
            ! Coarsen in x using 1-4-6-4-1 stencil:
            do ix = 0, cnx-1
                jx = 3 * ix
                jxm2 = mod(jx - 2 + fnx, fnx)
                jxm1 = mod(jx - 1 + fnx, fnx)
                jxp1 = mod(jx + 1, fnx)
                jxp2 = mod(jx + 2, fnx)

                do iy = 0, fny-1
                    do iz = 0, fnz
                        xdata(iz, iy, ix, :) = f116 * (fdata(iz, iy, jxm2, :) + &
                                                       fdata(iz, iy, jxp2, :))  &
                                             + f14  * (fdata(iz, iy, jxm1, :) + &
                                                       fdata(iz, iy, jxp1, :))  &
                                             + f38  *  fdata(iz, iy, jx, :)
                    enddo
                enddo
            enddo

            !------------------------------------------------------------------
            ! Coarsen in y using 1-4-6-4-1 stencil:
            do ix = 0, cnx-1
                do iy = 0, cny-1
                    jy = 3 * iy
                    jym2 = mod(jy - 2 + fny, fny)
                    jym1 = mod(jy - 1 + fny, fny)
                    jyp1 = mod(jy + 1, fny)
                    jyp2 = mod(jy + 2, fny)
                    do iz = 0, fnz
                        xydata(iz, iy, ix, :) = f116 * (xdata(iz, jym2, ix, :) + &
                                                        xdata(iz, jyp2, ix, :))  &
                                              + f14  * (xdata(iz, jym1, ix, :) + &
                                                        xdata(iz, jyp1, ix, :))  &
                                              + f38  *  xdata(iz, jy, ix, :)
                    enddo
                enddo
            enddo

            !------------------------------------------------------------------
            ! Coarsen in z using 1-4-6-4-1 stencil:
            ! Note: The boundaries are kept the same.
            do ix = 0, cnx-1
                do iy = 0, cny-1
                    do iz = 1, cnz - 1
                        jz = 3 * iz
                        cdata(iz, iy, ix, :) = f116 * (xydata(jz-2, iy, ix, :) + &
                                                       xydata(jz+2, iy, ix, :))  &
                                             + f14  * (xydata(jz-1, iy, ix, :) + &
                                                       xydata(jz+1, iy, ix, :))  &
                                             + f38  *  xydata(jz, iy, ix, :)
                    enddo
                    cdata(0,   iy, ix, :) = xydata(0,   iy, ix, :)
                    cdata(cnz, iy, ix, :) = xydata(fnz, iy, ix, :)
                enddo
            enddo
        end subroutine coarsen_three

        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer            :: i, stat
            character(len=512) :: arg
            logical            :: exists

            shrink = -1
            fname = ''
            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--filename') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    fname = trim(arg)
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
                    print *, 'This program coarsens all fields by a factor 2 or 3.'
                    print *, 'An EPIC field output must be provided to analyse.'
                    print *, 'Run code with "coarsen --filename [field file]" --shrink [coarsen shrink]'
                    stop
                endif
                i = i+1
            enddo

            if ((fname == '') .or. (shrink == -1) ) then
                print *, 'No file or step provided. Run code with "coarsen --help"'
                stop
            endif

            ! check if correct file is passed
            stat = index(trim(fname), '_fields.nc', back=.true.)
            if (stat == 0) then
                print *, "Error: No EPIC field output file provided."
                stop
            endif

            ! check if file exsits
            inquire(file=trim(fname), exist=exists)
            if (.not. exists) then
                print *, "Error: File '" // trim(fname) // "' does not exist."
                stop
            endif
        end subroutine parse_command_line

end program coarsening
