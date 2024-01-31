!   1. Determine local nearest interior/surface parcel to a small surface/interior parcel
!   2. Determine global nearest interior/surface parcel to a small surface/interior parcel
!   3. Send small interior/surface parcel to MPI rank with closest surface/interior parcel
!   4. Mix small parcel with its closest parcel
!   5. Send small parcel back to original MPI rank
module parcel_mixing
    use constants, only : one, zero, two
    use parcel_ops, only : get_delx
    use parcel_container, only : pc_type
    use parcels_mod, only : parcels, top_parcels, bot_parcels, ellipse_pc_type
    use parameters, only : dx, dxi, vcell, hli, lower, extent, acell    &
                         , nx, nz, vmin, max_num_parcels, amin          &
                         , max_num_surf_parcels
    use options, only : parcel
    use mpi_timer, only : start_timer, stop_timer
    use parcel_mpi, only : allocate_parcel_id_buffers       &
                         , deallocate_parcel_id_buffers
    use mpi_environment
    use mpi_utils, only : mpi_exit_on_error         &
                        , mpi_check_for_error       &
                        , mpi_check_for_message     &
                        , mpi_stop
    use parcel_nearest, only : nearest_type
    implicit none

    integer :: mixing_timer

    private

    !Used for searching for possible parcel mixer:
    logical, allocatable :: top_mixed(:), bot_mixed(:)

    integer, allocatable          :: isma(:)
    integer, allocatable          :: iclo(:)
    integer, allocatable          :: rclo(:)    ! MPI rank of closest parcel
    double precision, allocatable :: dclo(:)    ! distance to closest parcel


    type(nearest_type) :: near

    public :: mix_parcels, mixing_timer

    contains

        subroutine mix_parcels

            call start_timer(mixing_timer)

            ! 1. find small parcels in cells 0 and nz
            ! 2. find closest parcels
            ! 3. mix properties

            if (.not. allocated(top_mixed)) then
                allocate(top_mixed(max_num_surf_parcels))
                allocate(bot_mixed(max_num_surf_parcels))
            endif

            call interior2surface(0,  bot_parcels, bot_mixed)
            call interior2surface(nz, top_parcels, top_mixed)

            if (allocated(near%nppc)) then
                deallocate(near%nppc)
                deallocate(near%kc1)
                deallocate(near%kc2)
                deallocate(near%loca)
                deallocate(near%node)
            endif


            !------------------------------------------------------------------
            ! Mix small surface parcels with interior parcels:


            call find_interior_parcels

            call apply_mixing(bot_parcels, parcels, 0)

            call find_locally(bot_parcels, parcels, 0,  bot_mixed)

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(pcont, n_local_small, iclo, rclo, dclo)

            call actual_mixing(bot_parcels, parcels, isma, iclo, nmix)



            call find_locally(top_parcels, parcels, nz, top_mixed)
            call actual_mixing(top_parcels, parcels, isma, iclo, nmix)

            ! -----------------------------------------------------------------
            ! Mix small surface parcel with closest interior parcel:


            ! -----------------------------------------------------------------
            ! Final clean-up:

            deallocate(isma)
            deallocate(iclo)



            if (allocated(near%nppc)) then
                deallocate(near%nppc)
                deallocate(near%kc1)
                deallocate(near%kc2)
                deallocate(near%loca)
                deallocate(near%node)
            endif


            call stop_timer(mixing_timer)

        end subroutine mix_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine apply_mixing(source, dest, iz, is_mixed)
            class(pc_type),        intent(inout) :: source
            class(pc_type),        intent(inout) :: dest
            integer,               intent(in)    :: iz
            logical, allocatable,  intent(in)    :: is_mixed(:)
            integer                              :: n
            integer                              :: n_local_mix, n_global_mix, n_mix
            logical                              :: l_no_local_mix

            ! We must store the parcel index and the mix index *m*
            ! of each small parcel. We do not need to allocate the
            ! invalid buffer, therefore the second argument is .false.
            call allocate_parcel_id_buffers(source, 2, .false.)

            !------------------------------------------------------------------
            ! Initialise search:
            near%n_remote_small = 0

            ! -----------------------------------------------------------------
            ! Number of small surface parcels:
            n_local_mix = 0
            do n = 1, source%local_num
                if (source%is_small(n) .and. (.not. is_mixed(n))) then
                    n_local_mix = n_local_mix + 1

                    ! If a small parcel is in a boundary cell, a duplicate must
                    ! be sent to the neighbour rank. This call checks if the parcel
                    ! must be sent and fills the send buffers.
                    call near%locate_parcel_in_boundary_cell(source, n_local_mix, n)
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
                    "in MPI_Allreduce of parcel_mixing::surface2interior.")

            if (n_global_mix == 0) then
                call deallocate_parcel_id_buffers
                return
            endif

            !---------------------------------------------------------------------
            ! Communicate position of small parcels:
            ! Send position attribute, parcel index and the merge index of small parcels
            ! in boundary region to neighbours. We only need to send the position
            ! as this is the only attribute needed to figure out with whom a parcel
            ! might merge.
            call near%send_small_parcel_bndry_info(source)

            ! There are 4 cases:
            !   - n_local_mix = 0 and n_remote_small = 0 --> This rank has no small parcels.
            !   - n_local_mix > 0 and n_remote_small = 0 --> This rank has only local small parcels.
            !   - n_local_mix = 0 and n_remote_small > 0 --> This rank has only remote small parcels.
            !   - n_local_mix > 0 and n_remote_small > 0 --> This rank has local and remote small parcels.

            n_mix = n_local_mix + near%n_remote_small

            l_no_local_mix = (n_mix == 0)

            if (.not. l_no_local_mix) then

                allocate(isma(n_mix))
                allocate(iclo(n_mix))
                allocate(rclo(n_mix))
                allocate(dclo(n_mix))

                isma = -1
                iclo = -1
                rclo = -1
                dclo = product(extent)

                ! Fill isma array:
                m = 0
                do n = 1, source%local_num + near%n_remote_small
                    if (source%is_small(n) .and. (.not. is_mixed(n))) then
                        m = m + 1
                        isma(m) = n
                    endif
                enddo

                !---------------------------------------------------------------------
                ! Determine locally closest parcel:
                call find_locally(source, dest, iz, is_mixed)

            endif

            !---------------------------------------------------------------------
            ! Determine globally closest parcel:
            ! After this operation isma, iclo and rclo are properly set.
            call find_closest_parcel_globally(source, n_local_mix, iclo, rclo, dclo)

            if (.not. l_no_local_mix) then
                do n = 1, n_mix
                    if (iclo(n) == -1) then
                        ! ignore
                    endif
                enddo
            endif


            if (.not. l_no_local_mix) then
                !---------------------------------------------------------------------
                ! Mark all entries of isma, iclo and rclo above n_local_mix
                ! as invalid.
                do n = n_local_mix+1, size(iclo)
                    isma(n) = -1
                    iclo(n) = -1
                    rclo(n) = -1
                enddo

                isma(1:n_local_mix) = pack(isma(1:), isma(1:) /= -1)
                iclo(1:n_local_mix) = pack(iclo, iclo /= -1)
                rclo(1:n_local_mix) = pack(rclo, rclo /= -1)
            endif

            !---------------------------------------------------------------------
            ! We perform the actual merging locally. We must therefore send all
            ! necessary remote parcels to *this* MPI rank.
            ! Note: It cannot happen that the closest parcel is a small parcel
            !       on another MPI rank that is sent elsewhere.
            call gather_remote_parcels(pcont, n_local_small, n_invalid, rclo, iclo, isma, inva)


        end subroutine apply_mixing

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Find closest destination (dest) parcel to small source parcel:
        subroutine find_locally(source, dest, iz, is_mixed)
            class(pc_type),        intent(inout) :: source
            class(pc_type),        intent(inout) :: dest
            integer,               intent(in)    :: iz
            logical, allocatable,  intent(in)    :: is_mixed(:)
            integer                              :: n, ix, ij, m, ix0, ic, k, is
            double precision                     :: xs, ys, zs, delx, delz, dsq, dsqmin


            do m = 1, n_mix

                is = isma(m)

                ! position of surface parcel
                xs = source%position(1, is)
                ys = source%position(2, is)
                zs = source%get_z_position(is)

                ! grid cell this parcel is in:
                ix0 = nint(dxi(1) * (xs - box%halo_lower(1))) ! ranges from 0 to box%halo_size(1)
                iy0 = nint(dxi(2) * (ys - box%halo_lower(2))) ! ranges from 0 to box%halo_size(2)

                dsqmin = two * vcell
                ic = -1

                ! loop over all interior parcels in this grid cell
                ! and find closest parcel:
                do iy = max(0, iy0-1), min(box%size(2)+1, iy0)
                    do ix = max(0, ix0-1), min(box%size(1)+1, ix0)
                        ! Cell index:
                        ij = 1 + ix                                             &
                           + box%halo_size(1) * iy                              &
                           + box%halo_size(1) * box%halo_size(2) * min(iz, 1)

                        ! Search small parcels for closest other:
                        do k = near%kc1(ij), near%kc2(ij)
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
                isma(m) = is
                iclo(m) = ic
                rclo(m) = cart%rank
                dclo(m) = dsqmin

                ! If ic == -1, then no near interior parcel is found. We must ignore this mix.
                if (ic == -1) then
                    isma(m) = -1
                    rclo(m) = -1
                    dclo(m) = -1.0d0
                    n_mix = n_mix - 1
                endif
            enddo

            isma = pack(isma, isma /= -1)
            iclo = pack(iclo, iclo /= -1)
            rclo = pack(rclo, rclo /= -1)
            dclo = pack(dclo, dclo >= 0.0d0)

            !---------------------------------------------------------------------
            ! Update isma, iclo and rclo with indices of remote parcels:
            do m = 1, n_mix

                m = m + 1

                is = isma(n)
                ic = iclo(n)

                if ((is > source%local_num) .and. (ic > source%local_num)) then
                    ! A remote small parcel points to another remote small parcel. The remotes do not necessarily
                    ! need to be the same. Also, it can be a dual-link. As the same distance is evaluated on
                    ! the other two MPI ranks, *this* MPI rank must set the distance between the parcels to the
                    ! maximum value as otherwise the function "find_closest_parcel_globally" may think the
                    ! parcels belong to *this* MPI rank due to round-offs in the distance calculation. The
                    ! function "find_closest_parcel_globally" always sets "rclo" to the MPI source.
                    dclo(m) = huge(0.0d0) ! huge(x) returns the maximum value of this type
                else if (ic > source%local_num) then
                    ! A local small parcel points to a remote small parcel.
                    ! The index *ic* is larger than the local number of parcels
                    ! we must update *iclo(m)* and *rclo(m)* with the index of the parcel stored
                    ! on the other MPI rank.
                    iclo(m) = near%pidsma(ic)
                    rclo(m) = near%rsma(ic)
                endif
            enddo

            !------------------------------------------------------------------
            ! Sanity check: Indices of small parcels must be smaller equal to the
            ! number of local parcels:
            do m = 1, n_local_mix
                if (isma(m) > source%local_num) then
                    call mpi_exit_on_error(&
                        'in in parcel_mixing::find_locally: Small parcel index out of range.')
                endif
            enddo

        end subroutine find_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine find_interior_parcels
            integer :: iz, ij, n, ix, k

            if (.not. allocated(near%nppc)) then
                allocate(near%nppc(2 * nx))
                allocate(near%kc1(2 * nx))
                allocate(near%kc2(2 * nx))
                allocate(near%loca(max_num_parcels))
                allocate(near%node(max_num_parcels))
            endif

            !---------------------------------------------------------------------
            ! Initialise search:
            near%nppc = 0 !nppc(i) will contain the number of parcels in grid cell i

            ! Bin interior parcels in cells:
            do n = 1, parcels%local_num

                ! we only consider big parcels
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                iz = nint(dxi(2) * (parcels%position(2, n) - lower(2)))

                if ((iz == nz) .or. (iz == 0)) then
                    !
                    ! top or bottom surface
                    !

                    ix = mod(int(dxi(1) * (parcels%position(1, n) - lower(1))), nx)

                    ! Cell index of parcel:
                    ij = 1 + mod(nx + ix, nx) + nx * min(iz, 1)


                    ! Accumulate number of parcels in this grid cell:
                    near%nppc(ij) = near%nppc(ij) + 1

                    ! Store grid cell that this parcel is in:
                    near%loca(n) = ij
                endif
            enddo

            ! Find arrays kc1(i) & kc2(i) which indicate the parcels in grid cell i
            ! through n = node(k), for k = kc1(i),kc2(i):
            near%kc1(1) = 1
            do ij = 1, 2 * nx - 1
                near%kc1(ij+1) = near%kc1(ij) + near%nppc(ij)
            enddo

            near%kc2 = near%kc1 - 1
            do n = 1, parcels%local_num
                if (parcels%volume(n) < vmin) then
                    cycle
                endif

                iz = nint(dxi(2) * (parcels%position(2, n) - lower(2)))
                if ((iz == nz) .or. (iz == 0)) then
                    ij = near%loca(n)
                    k = near%kc2(ij) + 1
                    near%node(k) = n
                    near%kc2(ij) = k
                endif
            enddo

        end subroutine find_interior_parcels

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine actual_mixing(source, dest, isma, iclo, nmix)
            class(pc_type),        intent(inout) :: source
            class(pc_type),        intent(inout) :: dest
            integer,               intent(in)    :: isma(:)
            integer,               intent(in)    :: iclo(:)
            integer,               intent(in)    :: nmix
            integer                              :: l, m, n, ic, is
            double precision                     :: buoym(nmix), vortm(3, nmix), vm(nmix), vmix
            integer                              :: lclo(dest%local_num)

            lclo = zero

            !------------------------------------------------------------------
            ! Figure out weights:

            l = 0
            do m = 1, nmix

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

            do m = 1, nmix
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

end module parcel_mixing
