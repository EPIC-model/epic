program coarsening
    use constants, only : f14, f18, f16, f112
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
            integer                       :: v, dimids(4), varids(32), nc
            character(len=64)             :: name
            character(len=64)             :: fields(32)
            double precision, allocatable :: fdata(:, :, :, :)
            double precision, allocatable :: cdata(:, :, :, :)
            character(512)                :: fname
            integer                       :: coord_ids(3)  ! = (x, y, z)
            integer                       :: t_axis_id
            double precision              :: dx(3)
            double precision, allocatable :: t(:)

            call open_netcdf_file(trim(filename), NF90_NOWRITE, ncid, l_serial=.true.)

            call get_num_steps(ncid, n_steps)

            allocate(t(n_steps))

            ncerr = nf90_inquire(ncid, nDimensions, nVariables)

            ! allocate for all variables (excluding spatial and temporal variables)
            allocate(fdata(0:fnz, 0:fny-1, 0:fnx-1, nVariables - 4))
            allocate(cdata(0:cnz, 0:cny-1, 0:cnx-1, nVariables - 4))

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

            call define_netcdf_spatial_dimensions_3d(ncid=mcid,                &
                                                     ngps=(/cnx, cny, cnz+1/), &
                                                     dimids=dimids(1:3),       &
                                                     axids=coord_ids)

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

            call close_netcdf_file(mcid, l_serial=.true.)

            call open_netcdf_file(trim(fname), NF90_WRITE, mcid, l_serial=.true.)

            dx = extent / dble(ncells)
            call write_netcdf_axis_3d(mcid, dimids(1:3), lower, dx, &
                                      (/cnx, cny, cnz+1/))

            call read_netcdf_dataset(ncid, 't', t, start=(/1/), cnt=(/n_steps/))

            do step = 1, n_steps
                start = (/ 1,  1,  1, step /)

                nc = 1

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

                    call read_netcdf_dataset(ncid, trim(fields(v)), fdata(:, :, :, nc), &
                                             start=start, cnt=cnt)

                    nc = nc + 1
                enddo


                    if (shrink == 2) then
                        call coarsen_two(fdata, cdata, nVariables-4)
                    else if (shrink == 3) then
!                         print *, "Shrink by factor 3"
                        call coarsen_three(fdata, cdata, nVariables-4)
                    endif



                nc = 1
                do v = 1, nVariables


                    cnt = (/ cnx, cny, cnz+1, 1/)

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

                    call write_netcdf_dataset(mcid, v, cdata(:, :, :, nc), start=start, cnt=cnt, l_serial=.true.)

                    nc = nc + 1

                enddo
            enddo

            call close_netcdf_file(mcid, l_serial=.true.)
            call close_netcdf_file(ncid, l_serial=.true.)

            deallocate(t)

        end subroutine process_steps

        subroutine coarsen_two(fdata, cdata, nc)
            integer, intent(in) :: nc
            double precision, intent(in)  :: fdata(0:fnz, 0:fny-1, 0:fnx-1, nc)
            double precision, intent(out) :: cdata(0:cnz, 0:cny-1, 0:cnx-1, nc)
            double precision              :: xdata(0:fnz, 0:fny-1, 0:cnx-1, nc)
            double precision              :: xydata(0:fnz, 0:cny-1, 0:cnx-1, nc)
            integer                       :: ix, jx, iy, iz, jy, jz, jxp1, jyp1, jxm1, jym1


            do ix = 0, cnx-1
                jx = 2 * ix
                jxm1 = mod(jx - 1 + fnx, fnx)
                jxp1 = mod(jx + 1, fnx)

                do iy = 0, fny-1
                    do iz = 0, fnz
                        xdata(iz, iy, ix, :) = 1.0d0/4.0d0  * (fdata(iz, iy, jxm1, :) + fdata(iz, iy, jxp1, :)) &
                                          + 2.0d0/4.0d0 *  fdata(iz, iy, jx, :)
                    enddo
                enddo
            enddo

            do ix = 0, cnx-1
                do iy = 0, cny-1
                    jy = 2 * iy
                    jym1 = mod(jy - 1 + fny, fny)
                    jyp1 = mod(jy + 1, fny)
                    do iz = 0, fnz
                        xydata(iz, iy, ix, :) = 1.0d0/4.0d0  * (xdata(iz, jym1, ix, :) + xdata(iz, jyp1, ix, :)) &
                                           + 2.0d0/4.0d0 *  xdata(iz, jy, ix, :)
                    enddo
                enddo
            enddo

            do ix = 0, cnx-1
                do iy = 0, cny-1
                    do iz = 1, cnz - 1
                        jz = 2 * iz
                        cdata(iz, iy, ix, :) = 1.0d0/4.0d0  * (xydata(jz-1, iy, ix, :) + xydata(jz+1, iy, ix, :)) &
                                          + 2.0d0/4.0d0 *  xydata(jz, iy, ix, :)
                    enddo
                    cdata(0, iy, ix, :) = xydata(0, iy, ix, :)
                    cdata(cnz, iy, ix, :) = xydata(fnz, iy, ix, :)
                enddo
            enddo

!             do ix = 0, cnx-1
!                 jx = 2 * ix
!                 jxm1 = mod(jx - 1 + fnx, fnx)
!                 jxp1 = mod(jx + 1, fnx)
!                 do iy = 0, cny-1
!                     jy = 2 * iy
!                     jym1 = mod(jy - 1 + fny, fny)
!                     jyp1 = mod(jy + 1, fny)
!                     do iz = 1, cnz
!                         jz = 2 * iz - 1
!                         cdata(iz, iy, ix) = f112 * (fdata(jz, jy, jxm1) + fdata(jz, jy, jxp1)  &
!                                                   + fdata(jz, jym1, jx) + fdata(jz, jyp1, jx)  &
!                                                   + fdata(jz-1, jy, jx) + fdata(jz+1, jy, jx)) &
!                                           + f16 * fdata(jz, jy, jx)                           &
!                                           + f16 * fdata(jz, jy, jx)                           &
!                                           + f16 * fdata(jz, jy, jx)
!                     enddo
!                     cdata(0, iy, ix) = f18 * (fdata(0, jy, jxm1) + fdata(0, jy, jxp1)  &
!                                             + fdata(0, jym1, jx) + fdata(0, jyp1, jx)) &
!                                      + f14 * fdata(0, jy, jx)                          &
!                                      + f14 * fdata(0, jy, jx)
!
!                     cdata(cnz, iy, ix) = f18 * (fdata(fnz, jy, jxm1) + fdata(fnz, jy, jxp1)  &
!                                               + fdata(fnz, jym1, jx) + fdata(fnz, jyp1, jx)) &
!                                          + f14 * fdata(fnz, jy, jx)                          &
!                                          + f14 * fdata(fnz, jy, jx)
!                 enddo
!             enddo

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

            do ix = 0, cnx-1
                jx = 3 * ix
                jxm2 = mod(jx - 2 + fnx, fnx)
                jxm1 = mod(jx - 1 + fnx, fnx)
                jxp1 = mod(jx + 1, fnx)
                jxp2 = mod(jx + 2, fnx)

                do iy = 0, fny-1
                    do iz = 0, fnz
                        xdata(iz, iy, ix, :) = 1.0d0/16.d0  * (fdata(iz, iy, jxm2, :) + fdata(iz, iy, jxp2, :)) &
                                          + 4.0d0/16.d0  * (fdata(iz, iy, jxm1, :) + fdata(iz, iy, jxp1, :)) &
                                          + 6.0d0/16.0d0 *  fdata(iz, iy, jx, :)
                    enddo
                enddo
            enddo

            do ix = 0, cnx-1
                do iy = 0, cny-1
                    jy = 3 * iy
                    jym2 = mod(jy - 2 + fny, fny)
                    jym1 = mod(jy - 1 + fny, fny)
                    jyp1 = mod(jy + 1, fny)
                    jyp2 = mod(jy + 2, fny)
                    do iz = 0, fnz
                        xydata(iz, iy, ix, :) = 1.0d0/16.d0  * (xdata(iz, jym2, ix, :) + xdata(iz, jyp2, ix, :)) &
                                          + 4.0d0/16.d0  * (xdata(iz, jym1, ix, :) + xdata(iz, jyp1, ix, :)) &
                                          + 6.0d0/16.0d0 *  xdata(iz, jy, ix, :)
                    enddo
                enddo
            enddo

            do ix = 0, cnx-1
                do iy = 0, cny-1
                    do iz = 1, cnz - 1
                        jz = 3 * iz
                        cdata(iz, iy, ix, :) = 1.0d0/16.d0  * (xydata(jz-2, iy, ix, :) + xydata(jz+2, iy, ix, :)) &
                                          + 4.0d0/16.d0  * (xydata(jz-1, iy, ix, :) + xydata(jz+1, iy, ix, :)) &
                                          + 6.0d0/16.0d0 *  xydata(jz, iy, ix, :)
                    enddo
                    cdata(0, iy, ix, :) = xydata(0, iy, ix, :)
                    cdata(cnz, iy, ix, :) = xydata(fnz, iy, ix, :)
                enddo
            enddo

!
!
!
!
! !                 print *, ix, jxm2, jxm1, jx, jxp1, jxp2
!                 do iy = 0, cny-1
!                     jy = 3 * iy
!                     jym2 = mod(jy - 2 + fny, fny)
!                     jym1 = mod(jy - 1 + fny, fny)
!                     jyp1 = mod(jy + 1, fny)
!                     jyp2 = mod(jy + 2, fny)
!                     do iz = 1, cnz
!                         jz = 3 * iz - 1
! !                         cdata(iz, iy, ix) = 1.0d0 / 48.0d0 * (fdata(iz, jy, jxm2) + fdata(iz, jy, jxp2)   &
! !                                             +    fdata(iz, jym2, jx) + fdata(iz, jyp2, jx)   &
! !                                             +    fdata(iz-2, jy, jx) + fdata(iz+2, jy, jx))  &
! !                                           + 4.0d0 / 48.0d0 * (fdata(iz, jy, jxm1) + fdata(iz, jy, jxp1)   &
! !                                               +  fdata(iz, jym1, jx) + fdata(iz, jyp1, jx)   &
! !                                               +  fdata(iz-1, jy, jx) + fdata(iz+1, jy, jx))  &
! !                                           + 6.0d0 / 48.0d0 * fdata(iz, iy, jx)                            &
! !                                           + 6.0d0 / 48.0d0 * fdata(iz, iy, jx)                            &
! !                                           + 6.0d0 / 48.0d0 * fdata(iz, iy, jx)
!                     enddo
!                     cdata(0, iy, ix) = 1.0d0 / 32.0d0 * (fdata(0, jy, jxm2) + fdata(0, jy, jxp2)   &
!                                                     +     fdata(0, jym2, jx) + fdata(0, jyp2, jx))  &
!                                      + 4.0d0 / 32.0d0 * (fdata(0, jy, jxm1) + fdata(0, jy, jxp1)   &
!                                                      +    fdata(0, jym1, jx) + fdata(0, jyp1, jx))  &
!                                      + 6.0d0 / 32.0d0 * fdata(0, iy, jx)                           &
!                                      + 6.0d0 / 32.0d0 * fdata(0, iy, jx)
!
!
!
!                     cdata(cnz, iy, ix) = 1.0d0 / 32.0d0 * (fdata(fnz, jy, jxm2) + fdata(fnz, jy, jxp2)   &
!                                                      +    fdata(fnz, jym2, jx) + fdata(fnz, jyp2, jx))   &
!                                        + 4.0d0 / 32.0d0 * (fdata(fnz, jy, jxm1) + fdata(fnz, jy, jxp1)   &
!                                                      +    fdata(fnz, jym1, jx) + fdata(fnz, jyp1, jx))   &
!                                        + 6.0d0 / 32.0d0 * fdata(fnz, iy, jx)                           &
!                                        + 6.0d0 / 32.0d0 * fdata(fnz, iy, jx)
!                 enddo
!             enddo
!             stop


        end subroutine coarsen_three

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

end program coarsening
