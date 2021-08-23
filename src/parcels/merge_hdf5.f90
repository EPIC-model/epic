module merge_hdf5
    use parcel_container, only : parcels, n_parcels
    use hdf5
    use h5_utils
    use h5_writer
    use fields, only : get_index
    use merge_sort
    implicit none

    integer :: nw

    ! h5 file handles
    integer(hid_t)     :: h5file_id1, h5file_id2, h5file_id3
    character(len=512) :: h5fname1, h5fname2, h5fname3

    private :: h5file_id1, h5file_id2, h5fname1, h5fname2, h5fname3, nw

    contains

        subroutine create_h5_merger_files(basename, overwrite)
            character(*),   intent(in) :: basename
            logical,        intent(in) :: overwrite

            call create_h5_merger_file(h5file_id1, h5fname1, 'mergee', basename, overwrite)
            call create_h5_merger_file(h5file_id2, h5fname2, 'merger', basename, overwrite)
            call create_h5_merger_file(h5file_id3, h5fname3, 'neighbours', basename, overwrite)

        end subroutine create_h5_merger_files

        subroutine create_h5_merger_file(h5file_id, h5fname, name, basename, overwrite)
            integer(hid_t),     intent(inout) :: h5file_id
            character(len=512), intent(inout) :: h5fname
            character(*),       intent(in)    :: name
            character(*),       intent(in)    :: basename
            logical,            intent(in)    :: overwrite

            h5fname =  basename // '_' // name // '.hdf5'

            call create_h5_file(h5fname, overwrite, h5file_id)

            call write_h5_char_scalar_attrib(h5file_id, 'output_type', name)

            call write_h5_timestamp(h5file_id)
            call write_h5_box(h5file_id)

            call close_h5_file(h5file_id)

            nw = 0

        end subroutine create_h5_merger_file

        function get_group_merge_number(tag, i) result(name)
            character(*), intent(in) :: tag
            integer, intent(in)      :: i
            character(len=32)        :: name

            write(name, fmt='(I10.10)') i
            name = trim(tag) // '#' // trim(name)
        end function get_group_merge_number

        subroutine write_h5_parcels_in_cell(ibig, nmerge)
            integer, intent(in)           :: ibig(:)
            integer, intent(in)           :: nmerge
            integer(hid_t)                :: group, mgroup
            character(len=64)             :: tag
            character(:), allocatable     :: name
            logical                       :: created
            integer                       :: m, ib, n, nm, k, i, j, ii, jj, ilo, ihi, num
            integer                       :: ibig_sorted(nmerge), ind(nmerge)

            ibig_sorted = ibig(1:nmerge)

            ! inefficient to sort again, but it works
            call msort(ibig_sorted, ind)

            call open_h5_file(h5fname3, H5F_ACC_RDWR_F, h5file_id3)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id3, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id3, name, group)
            endif

            n = n_parcels + 1 ! make a gap to real parcels
            ib = -1
            nm = -1
            num = 0 ! merger number
            do m = 1, nmerge
                if (ib .ne. ibig_sorted(m)) then
                    ib = ibig_sorted(m)

                    call get_index(parcels%position(ib, :), i, j)
                    i = mod(i + nx, nx)
                    ilo = mod(i - 1 + nx, nx)
                    ihi = mod(i + 1 + nx, nx)

                    tag = get_group_merge_number('merger', num)

                    do k = 1, n_parcels
                        call get_index(parcels%position(k, :), ii, jj)
                        ii = mod(ii + nx, nx)

                        if ((j - 1 <= jj) .and. (jj < j + 2)) then
                            if ((ilo == ii) .or. (ii == i) .or. (ii == ihi)) then
                                    nm = nm + 1
                                    parcels%position(n + nm, :) = parcels%position(k, :)
                                    parcels%B(n + nm, :) = parcels%B(k, :)
                                    parcels%volume(n + nm) = parcels%volume(k)
                            endif
                        endif
                    enddo

                    call create_h5_group(group, trim(tag), mgroup)
                    ! write all involved parcels
                    call write_h5_dataset_2d(group, trim(tag), "position", &
                                             parcels%position(n:n + nm, :))

                    call write_h5_dataset_2d(group, trim(tag), "B", &
                                             parcels%B(n:n + nm, :))

                    call write_h5_dataset_1d(group, trim(tag), "volume", &
                                             parcels%volume(n:n + nm))

                    nm = -1
                    num = num + 1

                    call close_h5_group(mgroup)
                endif
            enddo

            call close_h5_group(group)

            call close_h5_file(h5file_id3)

        end subroutine write_h5_parcels_in_cell


        subroutine write_h5_mergees(isma, ibig, nmerge)
            integer, intent(in)           :: isma(0:)
            integer, intent(in)           :: ibig(:)
            integer, intent(in)           :: nmerge
            integer(hid_t)                :: group, mgroup
            character(len=64)             :: tag
            character(:), allocatable     :: name
            logical                       :: created
            integer                       :: m, n, l, ib, nm, is, num
            integer                       :: ibig_sorted(nmerge), ind(nmerge)

            ibig_sorted = ibig(1:nmerge)

            call msort(ibig_sorted, ind)

            call open_h5_file(h5fname1, H5F_ACC_RDWR_F, h5file_id1)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id1, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id1, name, group)
            endif


            n = n_parcels + 1 ! make a gap to real parcels
            nm = 0  ! number of involved parcels
            num = 0 ! merger number
            m = 1
            do l = 1, nmerge
                ib = ibig_sorted(l)

!                 print *, "ib = ", ib

                if (ib == 0) then
                    cycle
                endif

                ! big parcel
                parcels%position(n, :) = parcels%position(ib, :)
                parcels%B(n, :) = parcels%B(ib, :)
                parcels%volume(n) = parcels%volume(ib)

                nm = 1

                do m = l, nmerge
                    if (ib == ibig_sorted(m)) then
                        is = isma(ind(m))
                        parcels%position(n + nm, :) = parcels%position(is, :)
                        parcels%B(n + nm, :) = parcels%B(is, :)
                        parcels%volume(n + nm) = parcels%volume(is)
                        nm = nm + 1

!                         print *, "  is = ", is

                        ibig_sorted(m) = 0
                    endif
                enddo

!                 print *, "write", nm

                tag = get_group_merge_number('merger', num)

                call create_h5_group(group, trim(tag), mgroup)

                ! write all involved parcels
                call write_h5_dataset_2d(group, trim(tag), "position", &
                                          parcels%position(n:n + nm-1, :))

                call write_h5_dataset_2d(group, trim(tag), "B", &
                                         parcels%B(n:n + nm-1, :))

                call write_h5_dataset_1d(group, trim(tag), "volume", &
                                         parcels%volume(n:n+nm-1))

                call close_h5_group(mgroup)

                ! new big parcel
                num = num + 1

                ! reset number of involved parcels
                nm = 0
            enddo

            call write_h5_int_scalar_step_attrib(h5file_id1, nw, "nmergers", num)

            call close_h5_group(group)

            call close_h5_file(h5file_id1)

        end subroutine write_h5_mergees


        subroutine write_h5_mergers(ibig, nmerge)
            integer, intent(in)           :: ibig(:)
            integer, intent(in)           :: nmerge
            integer(hid_t)                :: group
            character(:), allocatable     :: name
            logical                       :: created
            integer                       :: m, ib, n, nm
            integer                       :: ibig_sorted(nmerge), ind(nmerge)

            ibig_sorted = ibig(1:nmerge)

            ! inefficient to sort again, but it works
            call msort(ibig_sorted, ind)

            call open_h5_file(h5fname2, H5F_ACC_RDWR_F, h5file_id2)

            name = trim(get_step_group_name(nw))

            call create_h5_group(h5file_id2, name, group, created)

            if (.not. created) then
                call open_h5_group(h5file_id2, name, group)
            endif

            n = n_parcels + 1 ! make a gap to real parcels
            ib = -1
            nm = -1
            do m = 1, nmerge
                if (ib .ne. ibig_sorted(m)) then
                    nm = nm + 1
                    ib = ibig_sorted(m)

                    parcels%position(n + nm, :) = parcels%position(ib, :)
                    parcels%B(n + nm, :) = parcels%B(ib, :)
                    parcels%volume(n + nm) = parcels%volume(ib)
                endif
            enddo


            ! write all involved parcels
            call write_h5_dataset_2d(h5file_id2, name, "position", &
                                     parcels%position(n:n + nm, :))

            call write_h5_dataset_2d(h5file_id2, name, "B", &
                                     parcels%B(n:n + nm, :))

            call write_h5_dataset_1d(h5file_id2, name, "volume", &
                                     parcels%volume(n:n + nm))

            call close_h5_group(group)

            call close_h5_file(h5file_id2)

!             print *, "done"
!             stop
            ! increase step number
            nw = nw + 1
        end subroutine write_h5_mergers
end module merge_hdf5
