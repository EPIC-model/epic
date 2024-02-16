!   1. Determine local nearest interior/surface parcel to a small surface/interior parcel
!   2. Determine global nearest interior/surface parcel to a small surface/interior parcel
!   3. Send small interior/surface parcel to MPI rank with closest surface/interior parcel
!   4. Mix small parcel with its closest parcel
!   5. Send small parcel back to original MPI rank
module parcel_mixing
    use constants, only : one, zero, two
    use parcel_ops, only : get_delx, get_dely
    use parcel_container, only : pc_type
    use parcels_mod, only : parcels, top_parcels, bot_parcels, ellipse_pc_type
    use parameters, only : dx, dxi, vcell, hli, lower, extent, acell    &
                         , nx, nz, vmin, max_num_parcels, amin          &
                         , max_num_surf_parcels
    use options, only : parcel
    use mpi_timer, only : start_timer, stop_timer
    use parcel_mpi, only : n_parcel_sends                   &
                         , allocate_parcel_id_buffers       &
                         , deallocate_parcel_id_buffers     &
                         , parcel_communicate
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error         &
                        , mpi_check_for_error       &
                        , mpi_check_for_message     &
                        , mpi_stop
    use mpi_layout, only : box, cart, neighbours    &
                         , get_neighbour_from_rank  &
                         , get_mpi_buffer           &
                         , allocate_mpi_buffers     &
                         , deallocate_mpi_buffers
    use parcel_nearest, only : near                             &
                             , locate_parcel_in_boundary_cell   &
                             , send_small_parcel_bndry_info     &
                             , update_remote_indices            &
                             , find_closest_parcel_globally
    implicit none

    integer :: mixing_timer

    private

    !Used for searching for possible parcel mixer:
    logical, allocatable :: top_mixed(:), bot_mixed(:), int_mixed(:)

    integer, allocatable          :: isma(:)
    integer, allocatable          :: inva(:)
    integer, allocatable          :: iclo(:)
    integer, allocatable          :: rclo(:)    ! MPI rank of closest parcel
    double precision, allocatable :: dclo(:)    ! distance to closest parcel

    integer :: n_parcel_recvs(8)

    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels

            call start_timer(mixing_timer)

            !------------------------------------------------------------------
            ! Initialise:

            if (.not. allocated(top_mixed)) then
                allocate(top_mixed(max_num_surf_parcels))
                allocate(bot_mixed(max_num_surf_parcels))
                allocate(int_mixed(max_num_parcels))
            endif


            ! We need to tag each parcel if it already
            ! mixed in order to avoid double-mixing:
            bot_mixed = .false.
            top_mixed = .false.
            int_mixed = .false.

            !------------------------------------------------------------------
            ! Mix small surface parcels with nearest interior parcels:

            call apply_mixing(bot_parcels, bot_mixed, parcels, int_mixed, 0)
            call apply_mixing(top_parcels, top_mixed, parcels, int_mixed, nz)

            ! -----------------------------------------------------------------
            ! Mix small interior parcels with nearest surface parcels:

            call apply_mixing(parcels, int_mixed, bot_parcels, bot_mixed, 0)
            call apply_mixing(parcels, int_mixed, top_parcels, top_mixed, nz)

            !------------------------------------------------------------------
            ! Final clean-up:

            if (allocated(top_mixed)) then
                deallocate(top_mixed)
                deallocate(bot_mixed)
                deallocate(int_mixed)
            endif

            call stop_timer(mixing_timer)

        end subroutine mix_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_mixing(source, l_sflag, dest, l_dflag, iz)
            class(pc_type), intent(inout) :: source
            logical,        intent(inout) :: l_sflag(:)
            class(pc_type), intent(inout) :: dest
            logical,        intent(inout) :: l_dflag(:)
            integer,        intent(in)    :: iz
            integer                       :: n, k
            integer                       :: n_local_mix, n_remote_mix
            integer                       :: n_global_mix, n_orig_parcels
            integer                       :: n_invalid, ic, is, m, j
            logical                       :: l_local_mix

            call near%alloc(max_num_parcels)

            near%nppc = 0

            ! We must store the parcel index and the merge index *m*
            ! of each small parcel. We do not need to allocate the
            ! invalid buffer, therefore the second argument is .false.
            call allocate_parcel_id_buffers(source, 2, .false.)


            ! -----------------------------------------------------------------
            ! Number of small parcels:
            n_local_mix = 0
            do n = 1, source%local_num
                if (source%is_small(n) .and. (.not. l_sflag(n))) then
                    k = nint(dxi(3) * (source%get_z_position(n) - box%lower(3)))
                    if (k == iz) then
                        n_local_mix = n_local_mix + 1

                        ! If a small parcel is in a boundary cell, a duplicate must
                        ! be sent to the neighbour rank. This call checks if the parcel
                        ! must be sent and fills the send buffers.
                        call locate_parcel_in_boundary_cell(source, n_local_mix, n)
                    endif
                endif
            enddo

            call MPI_Allreduce(n_local_mix,   &
                               n_global_mix,  &
                               1,             &
                               MPI_INTEGER,   &
                               MPI_SUM,       &
                               world%comm,    &
                               world%err)

            call mpi_check_for_error(world, &
                    "in MPI_Allreduce of parcel_mixing::apply_mixing.")


!             print *, "global mix:", n_global_mix
!             print *, "local mix: ", n_local_mix

            if (n_global_mix == 0) then
                call near%dealloc
                call deallocate_parcel_id_buffers
                return
            endif


            !------------------------------------------------------------------
            ! Assign destination (dest) parcels to grid cells:
            ! (only take parcels into account that have not yet been mixed)
            do n = 1, dest%local_num
                k = nint(dxi(3) * (dest%get_z_position(n) - box%lower(3)))
                if ((k == iz) .and. (.not. l_dflag(n))) then
                    call near%parcel_to_local_cell_index(dest, n)
                endif
            enddo

            call near%set_lbound

            do n = 1, dest%local_num
                k = nint(dxi(3) * (dest%get_z_position(n) - box%lower(3)))
                if ((k == iz) .and. (.not. l_dflag(n))) then
                    call near%set_ubound(n)
                endif
            enddo

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute, parcel index and the *isma* index of small parcels
            ! in boundary region to neighbours. We only need to send the position
            ! as this is the only attribute needed to figure out with whom a parcel
            ! might mix.
            call send_small_parcel_bndry_info(source, n_remote_mix)


!                     ! To print out the result enable the following lines:
!     do k = 0, world%size-1
!         if (k == world%rank) then
!             do n = 1, source%local_num + n_remote_mix
!                 j = 0
!                 if (n > source%local_num) then
!                     j = 1
!                 endif
! !                 is = isma(n)
! !                 ic = iclo(n)
!                 print *, j, source%position(1, n), source%position(2, n), int(source%buoyancy(n))
! !                                     source%position(1, ic), source%position(2, ic), int(source%buoyancy(ic))
!             enddo
!         endif
!         call MPI_Barrier(world%comm, world%err)
!     enddo
!     call mpi_stop



            l_local_mix = (n_local_mix + n_remote_mix > 0)

            if (l_local_mix) then

                !--------------------------------------------------------------
                ! Allocate working arrays:
                allocate(isma(0:n_local_mix+n_remote_mix))
!                 allocate(inva(0:n_local_mix+n_remote_mix))
                allocate(iclo(n_local_mix+n_remote_mix))
                allocate(rclo(n_local_mix+n_remote_mix))
                allocate(dclo(n_local_mix+n_remote_mix))

                isma = -1
                iclo = -1
                rclo = -1
!                 inva = -1
                dclo = product(extent)

!                 print *, "sizes:", source%local_num, n_remote_mix, n_local_mix
                j = 0
                do n = 1, source%local_num + n_remote_mix
                    if (source%is_small(n) .and. (.not. l_sflag(n))) then
                        k = nint(dxi(3) * (source%get_z_position(n) - box%lower(3)))
                        if (k == iz) then
                            j = j + 1
                            isma(j) = n
                        endif
                    endif
                enddo

                !--------------------------------------------------------------
                ! Determine locally closest parcel:
                call find_locally(source, dest, n_local_mix, n_remote_mix)
            endif

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(source, n_local_mix, iclo, rclo, dclo)

            if (l_local_mix) then
                !------------------------------------------------------------------
                ! Remove all parcels that did not find a near neighbour:
                isma(0) = 0 ! Ensure we do not delete the zero index (although not needed)
                n_local_mix = count(isma(1:) /= -1)
                isma(1:n_local_mix) = pack(isma(1:), isma(1:) /= -1)
                iclo(1:n_local_mix) = pack(iclo, iclo /= -1)
                rclo(1:n_local_mix) = pack(rclo, rclo /= -1)
            endif

            call deallocate_parcel_id_buffers

            if (n_local_mix == 0) then
                call near%dealloc
                deallocate(isma)
!                 deallocate(inva)
                deallocate(iclo)
                deallocate(rclo)
                deallocate(dclo)
                print *, "No parcel found a near neighbour."
                return
            endif

            !---------------------------------------------------------------------
            ! Mark all entries of isma, iclo and rclo above n_local_small
            ! as invalid.
            do n = n_local_mix+1, size(iclo)
                isma(n) = -1
                iclo(n) = -1
                rclo(n) = -1
            enddo


            !------------------------------------------------------------------
            ! Store original parcel number before we communicate:
            n_orig_parcels = source%local_num

            ! The number of small parcels that we receive is at most "n_remote_mix":
            allocate(inva(n_remote_mix))

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: It cannot happen that the closest parcel is a small parcel
            !       on another MPI rank that is sent elsewhere.
            ! The array *inva* contains all indices of small parcels that we sent to another
            ! MPI rank. We use *inva* to overwrite the parcels with their mixed values
            call gather_remote_parcels(source, n_local_mix, n_invalid)

            !------------------------------------------------------------------
            ! Apply the mixing and mark all parcels involved in the mixing
            ! process to avoid double-mixing
            if (n_local_mix > 0) then
                call actual_mixing(source, dest, n_local_mix)

                ! local mixing: applies to isma and iclo
                do m = 1, n_local_mix
                    is = isma(m)
                    ic = iclo(m)

                    ! Mark source parcel as 'mixed':
                    if (l_sflag(is)) then
                        call mpi_exit_on_error(&
                            'in apply_mixing: Small parcel was already mixed previously.')
                    endif
                    l_sflag(is) = .true.

                    ! Mark destination (dest) parcel as 'mixed':
                    if (l_dflag(ic)) then
                        call mpi_exit_on_error(&
                            'in apply_mixing: Close parcel was already mixed previously.')
                    endif
                    l_dflag(ic) = .true.
                enddo
            endif


            !------------------------------------------------------------------
            ! Receive small parcels that we sent earlier and mark as mixed:
            call scatter_remote_parcels(source, l_sflag, n_orig_parcels)

            deallocate(inva)

!             !------------------------------------------------------------------
!             ! The received parcels are all appended to the end of the container.
!             ! We must now overwrite their initial values with the new mixed values.
!             ! For this purpose we can use the *inva* array as it stores all invalid
!             ! parcel indices.
!             do m = 1, n_invalid
!                 is = inva(m)
!                 if (l_sflag(is)) then
!                     call mpi_exit_on_error(&
!                         'in apply_mixing: Small parcel was already mixed previously.')
!                 endif
!                 l_sflag(is) = .true.
!             enddo

!             ! Back fill the array with received parcels
!             call source%delete(inva, n_invalid)

            !------------------------------------------------------------------
            ! Deallocate working arrays:
            if (l_local_mix) then
                deallocate(isma)
!                 deallocate(inva)
                deallocate(iclo)
                deallocate(rclo)
                deallocate(dclo)
            endif

            call near%dealloc

        end subroutine apply_mixing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Find closest destination (dest) parcel to small source parcel:
        subroutine find_locally(source, dest, n_local_mix, n_remote_mix)
            class(pc_type),  intent(inout) :: source
            class(pc_type),  intent(inout) :: dest
            integer,         intent(in)    :: n_local_mix
            integer,         intent(in)    :: n_remote_mix
            integer                        :: n, ix, iy, ijk, m, ix0, iy0, iz0, ic, k, is
            double precision               :: xs, ys, zs, delx, dely, delz, dsq, dsqmin


            do m = 1, n_local_mix + n_remote_mix

                is = isma(m)

                ! position of surface parcel
                xs = source%position(1, is)
                ys = source%position(2, is)
                zs = source%get_z_position(is)

                ! nearest grid point:
                ix0 = nint(dxi(1) * (xs - box%halo_lower(1))) ! ranges from 0 to box%halo_size(1)
                iy0 = nint(dxi(2) * (ys - box%halo_lower(2))) ! ranges from 0 to box%halo_size(2)

                ! grid cells to search:
                iz0 = min(int(dxi(3) * (zs - box%lower(3))), nz-1)      ! ranges from 0 to nz-1

                dsqmin = two * vcell
                ic = 0

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do iy = max(0, iy0-1), min(box%size(2)+1, iy0)
                    do ix = max(0, ix0-1), min(box%size(1)+1, ix0)
                        ! Cell index:
                        ijk = 1 + ix                                        &
                                + box%halo_size(1) * iy                     &
                                + box%halo_size(1) * box%halo_size(2) * iz0

                        ! Search small parcels for closest other:
                        do k = near%kc1(ijk), near%kc2(ijk)
                            n = near%node(k)
                            delz = dest%get_z_position(n) - zs
                            if (delz * delz < dsqmin) then
                                delx = get_delx(dest%position(1, n), xs)
                                dely = get_dely(dest%position(2, n), ys)
                                ! Minimise dsqmin
                                dsq = delx ** 2 + dely ** 2 + delz ** 2
                                if (dsq < dsqmin) then
                                    dsqmin = dsq
                                    ic = n
                                endif
                            endif
                        enddo
                    enddo
                enddo

                ! Store the index of the parcel to be mixed with:
                isma(m)     = is
                iclo(m)     = ic
                rclo(m)     = cart%rank
                dclo(m)     = dsqmin

!                 print *, "m, is, ic", m, is, ic

                ! If ic == 0, then no near parcel is found. We must ignore this mix.
                if (ic == 0) then
                    isma(m) = -1
                    iclo(m) = -1
                    rclo(m) = -1
                    dclo(m) = huge(0.0d0)
                endif
            enddo

            !---------------------------------------------------------------------
            ! Make some checks:
            do m = 1, n_local_mix + n_remote_mix
                ic = iclo(m)

                if (ic > dest%local_num) then
                     ! A local small parcel points (is) to a remote close parcel (ic).
                     ! This cannot happen!
                     call mpi_exit_on_error(&
                        'in parcel_mixing::find_locally: Small parcel cannot point to remote close parcel!')
                endif
            enddo

            !------------------------------------------------------------------
            ! Sanity check: Indices of small parcels must be smaller equal to the
            ! number of local parcels:
            do m = 1, n_local_mix
                if (isma(m) > source%local_num) then
                    call mpi_exit_on_error(&
                        'in parcel_mixing::find_locally: Small parcel index out of range.')
                endif
            enddo


        end subroutine find_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine actual_mixing(source, dest, n_mix)
            class(pc_type),  intent(inout) :: source
            class(pc_type),  intent(inout) :: dest
            integer,         intent(in)    :: n_mix
            integer                        :: l, m, n, ic, is
            double precision               :: buoym(n_mix), vortm(3, n_mix), vm(n_mix), vmix
            integer                        :: lclo(dest%local_num)

            lclo = zero

            !------------------------------------------------------------------
            ! Figure out weights:

            l = 0
            do m = 1, n_mix

                ! Index of closest interior parcel
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    ! Start a new mixing parcel, indexed l:
                    l = l + 1
                    lclo(ic) = l

                    vm(l) = dest%volume(ic)


                    buoym(l) = vm(l) * dest%buoyancy(ic)

                    vortm(:, l) = vm(l) * dest%vorticity(:, ic)

                endif

                ! Sum up all the small surface parcel contributions:
                is = isma(m)
                n = lclo(ic)

                vm(n) = vm(n) + source%volume(is)

                buoym(n) = buoym(n) + source%volume(is) * source%buoyancy(is)

                vortm(:, n) = vortm(:, n) + source%volume(is) * source%vorticity(:, is)
            enddo

            !------------------------------------------------------------------
            ! Obtain the mixed values:
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmix = one / vm(m)

                buoym(m) = vmix * buoym(m)

                vortm(:, m) = vmix * vortm(:, m)
            enddo


            !------------------------------------------------------------------
            ! Apply the blended values:
            lclo = zero
            l = 0

            do m = 1, n_mix
                ic = iclo(m)

                if (lclo(ic) == 0) then
                    l = l + 1
                    lclo(ic) = l

                    dest%buoyancy(ic) = buoym(l)
                    dest%vorticity(:, ic) = vortm(:, l)
                endif

                is = isma(m)
                n = lclo(ic)

                source%buoyancy(is) = buoym(l)
                source%vorticity(:, is) = vortm(:, l)
            enddo

        end subroutine actual_mixing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine gather_remote_parcels(pcont, n_local_small, n_invalid)
            class(pc_type),              intent(inout) :: pcont
            integer,                     intent(inout) :: n_local_small
            integer,                     intent(inout) :: n_invalid
            double precision, dimension(:), pointer    :: send_buf
            double precision, allocatable              :: recv_buf(:)
            type(MPI_Request)                          :: requests(8)
            type(MPI_Status)                           :: recv_status, send_statuses(8)
            integer                                    :: recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(pcont%mix_attr_num)
            integer                                    :: m, rc, ic, is, n, i, j, k, iv
            integer                                    :: n_entries
            integer                                    :: n_registered(8)

            !------------------------------------------------------------------
            ! Figure out how many parcels we send and allocate buffers:
            n_parcel_sends = 0

            do m = 1, n_local_small
                rc = rclo(m)
                if (cart%rank /= rc) then
                    n = get_neighbour_from_rank(rc)
                    n_parcel_sends(n) = n_parcel_sends(n) + 1
                endif
            enddo

            n_invalid = sum(n_parcel_sends)

            ! We must send all parcel attributes used for mixing (mix_attr_num), the local
            ! parcel index (is) plus the index of the close parcel ic (1)
            n_entries = pcont%mix_attr_num + 2
            n_registered = n_parcel_sends * n_entries
            call allocate_mpi_buffers(n_registered)

            !------------------------------------------------------------------
            ! Fill all buffers:
            n_registered = 0
            iv = 0
            do m = 1, n_local_small
                is = isma(m)
                ic = iclo(m)
                rc = rclo(m)
                if (cart%rank /= rc) then
                    ! The closest parcel to this small parcel *is*
                    ! is on another MPI rank. We must send this parcel
                    ! to that rank.

                    n = get_neighbour_from_rank(rc)

                    call get_mpi_buffer(n, send_buf)

                    n_registered(n) = n_registered(n) + 1

                    call pcont%mixing_serialize(is, buffer)
                    j = n_entries * n_registered(n)
                    i = j - n_entries + 1
                    send_buf(i:j-2) = buffer
                    send_buf(j-1) = dble(ic)
                    send_buf(j)   = dble(is)

                    ! Mark this as invalid
                    isma(m) = -1
                    iclo(m) = -1
                    rclo(m) = -1

                    iv = iv + 1
                endif
            enddo

            if (iv /= n_invalid) then
                call mpi_exit_on_error(&
                    "in parcel_mixing::gather_remote_parcels: Numbers of invalid parcels disagree.")
            endif

            if (n_local_small > 0) then
                ! We must now remove all invalid entries in isma and
                ! iclo and also update the value of n_local_small:
                iv = n_local_small
                n_local_small = count(isma(1:iv) /= -1)
                isma(1:n_local_small) = pack(isma(1:iv), isma(1:iv) /= -1)
                iclo(1:n_local_small) = pack(iclo(1:iv), iclo(1:iv) /= -1)
                rclo(1:n_local_small) = pack(rclo(1:iv), rclo(1:iv) /= -1)
            endif

            !------------------------------------------------------------------
            ! Communicate parcels:
            do n = 1, 8
                call get_mpi_buffer(n, send_buf)

                send_size = n_parcel_sends(n) * n_entries

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of parcel_mixing::gather_remote_parcels.")
            enddo

            iv = 0
            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size,               &
                                           cart)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              recv_status,              &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Recv of parcel_mixing::gather_remote_parcels.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_mixing::gather_remote_parcels: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                n_parcel_recvs(n) = recv_count

                if (recv_count > 0) then
                    ! Set the current value of the *m* index to
                    ! its current maximum value (excluding small
                    ! parcel received earlier from other MPI ranks.
                    m = n_local_small
                    do k = 1, recv_count
                        ! We receive a small parcel;
                        ! append it to the container with index
                        ! "pcont%local_num+1"
                        pcont%local_num = pcont%local_num + 1
                        j = n_entries * k
                        i = j - n_entries + 1
                        buffer = recv_buf(i:j-2)
                        call pcont%mixing_deserialize(pcont%local_num, buffer)

                        ! Add the small parcel to isma, its index to inva
                        ! and its closest parcel index to iclo
                        m = m + 1
                        iclo(m) = nint(recv_buf(j-1))
                        isma(m) = pcont%local_num
                        iv = iv + 1
                        inva(iv) = nint(recv_buf(j))
                    enddo
                    ! The last value of *m* is our new number of
                    ! small parcels. However, this still includes
                    ! invalid entries in isma and iclo that we will
                    ! remove next.
                    n_local_small = m
                endif

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of parcel_mixing::gather_remote_parcels.")

            call deallocate_mpi_buffers

#ifndef NDEBUG
            n = pcont%local_num - n_invalid
            call mpi_blocking_reduce(n, MPI_SUM, world)
            if ((world%rank == world%root) .and. (.not. n == pcont%total_num)) then
                call mpi_exit_on_error(&
                    "in parcel_mixing::gather_remote_parcels: We lost parcels.")
            endif
#endif
        end subroutine gather_remote_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Receive small parcels with remote *iclo* parcel from *this* remote MPI rank.
        ! The small parcels *this* MPI rank received earlier are all at the end of the
        ! container. The parcel indices are stored in *inva*.
        subroutine scatter_remote_parcels(pcont, l_flag, n_orig_parcels)
            class(pc_type),              intent(inout) :: pcont
            logical,                     intent(inout) :: l_flag(:)
            integer,                     intent(in)    :: n_orig_parcels
            double precision, dimension(:), pointer    :: send_buf
            double precision, allocatable              :: recv_buf(:)
            type(MPI_Request)                          :: requests(8)
            type(MPI_Status)                           :: recv_status, send_statuses(8)
            integer                                    :: recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(pcont%mix_attr_num)
            integer                                    :: m, n, i, j, k, iv
            integer                                    :: n_entries
            integer                                    :: n_registered(8)

            ! We must send all parcel attributes used for mixing (mix_attr_num) plus its
            ! original local parcel index, stored in *inva*.
            ! "n_parcel_recvs" is the number of parcels we received earlier, we must send
            ! these now back to their original MPI rank.
            n_entries = pcont%mix_attr_num + 1
            n_registered = n_parcel_recvs * n_entries
            call allocate_mpi_buffers(n_registered)

            !------------------------------------------------------------------
            ! Communicate parcels:

            ! start index of remote parcels
            k = 0
            do n = 1, 8
                call get_mpi_buffer(n, send_buf)

                ! Fill send buffer:
                do iv = 1, n_parcel_recvs(n)
                    m = k + iv + n_orig_parcels
                    call pcont%mixing_serialize(m, buffer)
                    j = n_entries * iv
                    i = j - n_entries + 1
                    send_buf(i:j-1) = buffer
                    send_buf(j)     = dble(inva(k+iv))
                enddo
                k = k + n_parcel_recvs(n)

                send_size = n_parcel_recvs(n) * n_entries

                call MPI_Isend(send_buf(1:send_size),   &
                               send_size,               &
                               MPI_DOUBLE_PRECISION,    &
                               neighbours(n)%rank,      &
                               SEND_NEIGHBOUR_TAG(n),   &
                               cart%comm,               &
                               requests(n),             &
                               cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Isend of parcel_mixing::scatter_remote_parcels.")
            enddo

            ! We need to reduce the parcel number count again:
            pcont%local_num = pcont%local_num - sum(n_parcel_recvs)

            do n = 1, 8
                ! check for incoming messages
                call mpi_check_for_message(neighbours(n)%rank,      &
                                           RECV_NEIGHBOUR_TAG(n),   &
                                           recv_size,               &
                                           cart)

                allocate(recv_buf(recv_size))

                call MPI_Recv(recv_buf(1:recv_size),    &
                              recv_size,                &
                              MPI_DOUBLE_PRECISION,     &
                              neighbours(n)%rank,       &
                              RECV_NEIGHBOUR_TAG(n),    &
                              cart%comm,                &
                              recv_status,              &
                              cart%err)

                call mpi_check_for_error(cart, &
                    "in MPI_Recv of parcel_mixing::scatter_remote_parcels.")

                if (mod(recv_size, n_entries) /= 0) then
                    call mpi_exit_on_error(&
                        "parcel_mixing::scatter_remote_parcels: Receiving wrong count.")
                endif

                recv_count = recv_size / n_entries

                if (recv_count > 0) then
                    do k = 1, recv_count
                        ! We receive a small mixed parcel;
                        ! replace its original values at index "iv"
                        ! with the new mixed values
                        j = n_entries * k
                        i = j - n_entries + 1
                        buffer = recv_buf(i:j-1)
                        iv = nint(recv_buf(j))
                        call pcont%mixing_deserialize(iv, buffer)

                        if (l_flag(iv)) then
                            call mpi_exit_on_error(&
                            'in parcel_mixing::scatter_remote_parcels: Small parcel was already mixed previously.')
                        endif
                        l_flag(iv) = .true.
                        print *, "receive parcel", iv
                    enddo
                endif

                deallocate(recv_buf)
            enddo

            call MPI_Waitall(8,                 &
                            requests,           &
                            send_statuses,      &
                            cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Waitall of parcel_mixing::scatter_remote_parcels.")

            call deallocate_mpi_buffers

#ifndef NDEBUG
            n = pcont%local_num
            call mpi_blocking_reduce(n, MPI_SUM, world)
            if ((world%rank == world%root) .and. (.not. n == pcont%total_num)) then
                call mpi_exit_on_error(&
                    "in parcel_mixing::scatter_remote_parcels: We lost parcels.")
            endif
#endif

        end subroutine scatter_remote_parcels

!         !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!         ! Receive small parcels with remote *iclo* parcel from *this* remote MPI rank.
!         ! The small parcels *this* MPI rank received earlier are all at the end of the
!         ! container. The parcel indices are stored in *inva*.
!         subroutine scatter_remote_parcels(source, n_orig_parcels)
!             class(pc_type), intent(inout) :: source
!             integer,        intent(in)    :: n_orig_parcels
!             integer, allocatable          :: pindex(:)
!             integer                       :: n, num
!
!             num = source%local_num - n_orig_parcels
!
!             if (num < 0) then
!                 call mpi_exit_on_error(&
!                     "in parcel_mixing::scatter_remote_parcels: We have fewer parcels as before.")
!             endif
!
!             allocate(pindex(num))
!
!             ! All parcels that *this* MPI rank must send are stored at indices
!             ! > n_orig_parcels:
!             do n = 1, num
!                 pindex(n) = n_orig_parcels + n
!             enddo
!
!             ! This call also deletes the parcels *this* MPI rank sends.
!             call parcel_communicate(source, pindex)
!
!             deallocate(pindex)
!
!         end subroutine scatter_remote_parcels

end module parcel_mixing
