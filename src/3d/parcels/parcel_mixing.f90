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
    use parcel_nearest, only : handle_periodic_edge_parcels     &
                             , locate_parcel_in_boundary_cell   &
                             , send_small_parcel_bndry_info     &
                             , update_remote_indices            &
                             , find_closest_parcel_globally
    implicit none

    integer :: mixing_timer

    private

    !Used for searching for possible parcel merger:
    integer, allocatable :: nppc(:), kc1(:), kc2(:)
    integer, allocatable :: loca(:)
    integer, allocatable :: node(:)

    !Used for searching for possible parcel mixer:
    logical, allocatable :: top_mixed(:), bot_mixed(:), int_mixed(:)

    integer, allocatable          :: isma(:)
    integer, allocatable          :: inva(:)
    integer, allocatable          :: iclo(:)
    integer, allocatable          :: rclo(:)    ! MPI rank of closest parcel
    double precision, allocatable :: dclo(:)    ! distance to closest parcel


    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels
            integer :: nelem

            call start_timer(mixing_timer)

            !------------------------------------------------------------------
            ! Initialise:

            if (.not. allocated(top_mixed)) then
                nelem = box%halo_size(1) * box%halo_size(2)
                allocate(nppc(nelem))
                allocate(kc1(nelem))
                allocate(kc2(nelem))
                allocate(loca(max_num_parcels))
                allocate(node(max_num_parcels))
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

            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
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
            integer                       :: n_mix, n_invalid, ic, is, m
            logical                       :: l_local_mix

            ! -----------------------------------------------------------------
            ! Number of small parcels:
            n_local_mix = 0
            do n = 1, source%local_num
                if (source%is_small(n) .and. (.not. l_sflag(n))) then
                    k = nint(dxi(3) * (source%get_z_position(n) - box%lower(3)))
                    if (k == iz) then
                        n_local_mix = n_local_mix + 1
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

            if (n_global_mix == 0) then
                return
            endif

            !------------------------------------------------------------------
            ! Assign destination (dest) parcels to grid cells:
            ! (only take parcels into account that have not yet been mixed)
            do n = 1, dest%local_num
                k = nint(dxi(3) * (dest%get_z_position(n) - box%lower(3)))
                if ((k == iz) .and. (.not. l_dflag(n))) then
                    call local_cell_index(dest, n)
                endif
            enddo

            ! Sets n_remote_mix
            call exchange_bndry_info(source, n_local_mix, n_remote_mix)

            l_local_mix = (n_local_mix + n_remote_mix > 0)

            if (l_local_mix) then

                !--------------------------------------------------------------
                ! Allocate working arrays:
                allocate(isma(0:n_local_mix+n_remote_mix))
                allocate(inva(0:n_local_mix+n_remote_mix))
                allocate(iclo(n_local_mix+n_remote_mix))
                allocate(rclo(n_local_mix+n_remote_mix))
                allocate(dclo(n_local_mix+n_remote_mix))

                isma = -1
                iclo = -1
                rclo = -1
                inva = -1
                dclo = product(extent)

                !--------------------------------------------------------------
                ! Determine locally closest parcel:
                call find_locally(source, dest, n_local_mix, n_remote_mix)
            endif

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(source, n_local_mix, iclo, rclo, dclo)

            !------------------------------------------------------------------
            ! Remove all parcels that did not find a near neighbour:
            isma(0) = 0 ! Ensure we do not delete the zero index (although not needed)
            isma = pack(isma, isma /= -1)
            iclo = pack(iclo, iclo /= -1)
            rclo = pack(rclo, rclo /= -1)

!             ! Count how many parcels *this* MPI rank sends to the remote:
!             ! Because rclo is an allocatable array, the call to pack shrinks
!             ! the array size; i.e. it is safe to count like this; this would not
!             ! work with an array not marked as *allocatable*.
!             n_remote = count(rclo /= cart%rank)

            !------------------------------------------------------------------
            ! Store original parcel number before we communicate:
            n_orig_parcels = source%local_num

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: It cannot happen that the closest parcel is a small parcel
            !       on another MPI rank that is sent elsewhere.
            ! The array *inva* contains all indices of small parcels that we sent to another
            ! MPI rank. We use *inva* to overwrite the parcels with their mixed values
            call gather_remote_parcels(source, n_local_mix, n_invalid, rclo, iclo, isma, inva)

!             !------------------------------------------------------------------
!             ! We cannot use *inva* as it is filled in "gather_remote_parcels"
!             ! directly for parcel deletion as it adds all small parcels for deletion
!             ! because this algorithm is used in the parcel merging where finally all small
!             ! parcel must be deleted. Here, we must ensure we only delete the
!             ! parcels that we sent. All parcels at the end of *inva* are received. After the
!             ! mixing we must send all parcels with index > n_orig_parcels back to the original
!             ! MPI rank. "n_invalid" contains the number of parcels *this* MPI rank sent.


            !------------------------------------------------------------------
            ! Apply the mixing:
            if (l_local_mix) then
                n_mix = n_local_mix !+ n_remote_mix
                call actual_mixing(source, dest, n_mix)
            endif

            !------------------------------------------------------------------
            ! Mark all parcels involved in the mixing process to
            ! avoid double-mixing:
            if (l_local_mix) then
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
            ! Receive small parcels that we sent earlier:
            call scatter_remote_parcels(source, n_orig_parcels)

            !------------------------------------------------------------------
            ! The received parcels are all appended to the end of the container.
            ! We must now overwrite their initial values with the new mixed values.
            ! For this purpose we can use the *inva* array as it stores all invalid
            ! parcel indices.
            do m = 1, n_invalid
                is = inva(m)
                if (l_sflag(is)) then
                    call mpi_exit_on_error(&
                        'in apply_mixing: Small parcel was already mixed previously.')
                endif
                l_sflag(is) = .true.
            enddo

            ! Back fill the array with received parcels
            call source%delete(inva, n_invalid)

            !------------------------------------------------------------------
            ! Deallocate working arrays:
            if (l_local_mix) then
                deallocate(isma)
                deallocate(inva)
                deallocate(iclo)
                deallocate(rclo)
                deallocate(dclo)
            endif

        end subroutine apply_mixing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Assign a parcel to (ix, iy). The vertical direction is not needed.
        !@pre Assumes a parcel is in the local domain including halo cells
        !      (in x and y).
        subroutine local_cell_index(pcont, n)
            class(pc_type), intent(inout) :: pcont
            integer,        intent(in)    :: n
            integer                       :: ix, iy, ij

            call handle_periodic_edge_parcels(pcont%position(:, n))

            ix = int(dxi(1) * (pcont%position(1, n) - box%halo_lower(1)))
            iy = int(dxi(2) * (pcont%position(2, n) - box%halo_lower(2)))

            ! Cell index of parcel:
            !   This runs from 1 to halo_ncell where
            !   halo_ncell includes halo cells
            ij = 1 + ix + box%halo_size(1) * iy

#ifndef NEDBUG
            if (ij < 1) then
               call mpi_exit_on_error(&
                        'in local_cell_index: Parcel not in local domain.')
            endif
#endif

            ! Accumulate number of parcels in this grid cell:
            nppc(ij) = nppc(ij) + 1

            ! Store grid cell that this parcel is in:
            loca(n) = ij

        end subroutine local_cell_index

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! If a small parcel is in a boundary cell, a duplicate must
        ! be sent to the neighbour rank.
        subroutine exchange_bndry_info(source, n_local_mix, n_remote_mix)
            class(pc_type), intent(inout) :: source
            integer,        intent(in)    :: n_local_mix
            integer,        intent(inout) :: n_remote_mix
            integer                       :: m, is

            ! We must store the parcel index and the merge index *m*
            ! of each small parcel. We do not need to allocate the
            ! invalid buffer, therefore the second argument is .false.
            call allocate_parcel_id_buffers(source, 2, .false.)

            !------------------------------------------------------------------
            ! Collect info and fill send buffers
            do m = 1, n_local_mix
                is = isma(m)
                call locate_parcel_in_boundary_cell(source, m, is)
            enddo

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute, parcel index and the *isma* index of small parcels
            ! in boundary region to neighbours. We only need to send the position
            ! as this is the only attribute needed to figure out with whom a parcel
            ! might mix.
            call send_small_parcel_bndry_info(source, n_remote_mix)

            !------------------------------------------------------------------
            ! Clean-up:
            call deallocate_parcel_id_buffers

        end subroutine exchange_bndry_info

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Find closest destination (dest) parcel to small source parcel:
        subroutine find_locally(source, dest, n_local_mix, n_remote_mix)
            class(pc_type),  intent(inout) :: source
            class(pc_type),  intent(inout) :: dest
            integer,         intent(in)    :: n_local_mix
            integer,         intent(in)    :: n_remote_mix
            integer                        :: n, ix, iy, ij, m, ix0, iy0, ic, k, is
            double precision               :: xs, ys, zs, delx, dely, delz, dsq, dsqmin


            do m = 1, n_local_mix + n_remote_mix

                is = isma(m)

                ! position of surface parcel
                xs = source%position(1, is)
                ys = source%position(2, is)
                zs = source%get_z_position(is)

                ! grid cell this parcel is in:
                ix0 = nint(dxi(1) * (xs - box%halo_lower(1))) ! ranges from 0 to box%halo_size(1)
                iy0 = nint(dxi(2) * (ys - box%halo_lower(2))) ! ranges from 0 to box%halo_size(2)

                dsqmin = two * vcell
                ic = 0

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do iy = max(0, iy0-1), min(box%size(2)+1, iy0)
                    do ix = max(0, ix0-1), min(box%size(1)+1, ix0)
                        ! Cell index:
                        ij = 1 + ix + box%halo_size(1) * iy

                        ! Search small parcels for closest other:
                        do k = kc1(ij), kc2(ij)
                            n = node(k)
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

                ! If ic == -1, then no near parcel is found. We must ignore this mix.
                if (ic == 0) then
                    isma(m) = -1
                    rclo(m) = -1
                    dclo(m) = huge(0.0d0)
                endif
            enddo

            !---------------------------------------------------------------------
            ! Update isma, iclo and rclo with indices of remote parcels:
            call update_remote_indices(source, n_local_mix, n_remote_mix, isma, iclo, rclo, dclo)

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

        subroutine gather_remote_parcels(pcont, n_local_small, n_invalid, rclo, iclo, isma, inva)
            class(pc_type),              intent(inout) :: pcont
            integer,                     intent(inout) :: n_local_small
            integer,                     intent(inout) :: n_invalid
            integer,                     intent(inout) :: rclo(:)       ! MPI rank of closest parcel
            integer,                     intent(inout) :: iclo(:)
            integer,                     intent(inout) :: isma(0:)
            integer,                     intent(inout) :: inva(0:)
            double precision, dimension(:), pointer    :: send_buf
            double precision, allocatable              :: recv_buf(:)
            type(MPI_Request)                          :: requests(8)
            type(MPI_Status)                           :: recv_status, send_statuses(8)
            integer                                    :: recv_count
            integer                                    :: recv_size, send_size
            double precision                           :: buffer(pcont%attr_num)
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

            ! We must send all parcel attributes (attr_num) plus
            ! the index of the close parcel ic (1)
            n_entries = pcont%attr_num + 1
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

                    call pcont%serialize(is, buffer)
                    j = n_entries * n_registered(n)
                    i = j - n_entries + 1
                    send_buf(i:j-1) = buffer
                    send_buf(j) = dble(ic)

                    ! Mark this as invalid
                    isma(m) = -1
                    iclo(m) = -1
                    rclo(m) = -1

                    ! Copy index of this parcel, we will receive this parcel back
                    iv = iv + 1
                    inva(iv) = is
                endif
            enddo

            if (iv /= n_invalid) then
                call mpi_exit_on_error(&
                    "in parcel_mixing::gather_remote_parcels: Numbers of invalid parcels disagree.")
            endif

            if (n_local_small > 0) then
                ! We must now remove all invalid entries in isma and
                ! iclo and also update the value of n_local_small:
                n_local_small = count(isma(1:) /= -1)
                isma(1:n_local_small) = pack(isma(1:), isma(1:) /= -1)
                iclo(1:n_local_small) = pack(iclo, iclo /= -1)
                rclo(1:n_local_small) = pack(rclo, rclo /= -1)
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
                        buffer = recv_buf(i:j-1)
                        call pcont%deserialize(pcont%local_num, buffer)

                        ! Add the small parcel to isma and
                        ! its closest parcel to iclo
                        m = m + 1
                        iclo(m) = nint(recv_buf(j))
                        isma(m) = pcont%local_num
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
        ! container.
        subroutine scatter_remote_parcels(source, n_orig_parcels)
            class(pc_type), intent(inout) :: source
            integer,        intent(in)    :: n_orig_parcels
            integer, allocatable          :: pindex(:)
            integer                       :: n, num

            num = source%local_num - n_orig_parcels

            if (num < 0) then
                call mpi_exit_on_error(&
                    "in parcel_mixing::scatter_remote_parcels: We have fewer parcels as before.")
            endif

            allocate(pindex(num))

            ! All parcels that *this* MPI rank must send are stored at indices
            ! > n_orig_parcels:
            do n = 1, num
                pindex(n) = n_orig_parcels + n
            enddo

            ! This call also deletes the parcels *this* MPI rank sends.
            call parcel_communicate(source, pindex)

            deallocate(pindex)

        end subroutine scatter_remote_parcels

end module parcel_mixing
