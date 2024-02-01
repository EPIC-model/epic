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
    use parcel_nearest, only : handle_periodic_edge_parcels
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
    integer, allocatable          :: iclo(:)
    integer, allocatable          :: rclo(:)    ! MPI rank of closest parcel
    double precision, allocatable :: dclo(:)    ! distance to closest parcel


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
                allocate(int_mixed(max_num_parcels))
            endif


            call interior2surface(0,  bot_parcels, bot_mixed)
            call interior2surface(nz, top_parcels, top_mixed)

            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
            endif


            !------------------------------------------------------------------
            ! Mix small surface parcels with interior parcels:

            ! We need to tag each small surface parcel if it already
            ! mixed with an interior parcel in order to avoid double-mixing:
            bot_mixed = .false.
            top_mixed = .false.
            int_mixed = .false.

            call apply_mixing(bot_parcels, parcels, 0)

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



            if (allocated(nppc)) then
                deallocate(nppc)
                deallocate(kc1)
                deallocate(kc2)
                deallocate(loca)
                deallocate(node)
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
            integer                       :: n_local_mix, n_global_mix, n_mix
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

            do n = 1, dest%local_num
                k = nint(dxi(3) * (dest%get_z_position(n) - box%lower(3)))
                if (k == iz) then
                    call local_cell_index(dest, n)
                endif
            enddo

            call find_locally(source, dest, l_dflag)

            call find_globally

            ! Send small parcels to MPI rank owning the close parcel
            call send_small_parcels

            call actual_mixing(source, dest, n_mix)

            ! Receive small parcels that we sent earlier
            call receive_small_parcels

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

        ! Find closest destination (dest) parcel to small source parcel:
        subroutine find_locally(source, dest, l_dflag)
            class(pc_type),  intent(inout) :: source
            class(pc_type),  intent(inout) :: dest
            logical,         intent(inout) :: l_dflag(:)
            integer                        :: n, ix, ij, m, ix0, ic, k, is
            double precision               :: xs, ys, zs, delx, delz, dsq, dsqmin


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
                    n_mix = n_mix - 1
                else
                    l_dflag(ic) = .true.
                endif
            enddo

!             isma = pack(isma, isma /= -1)
!             iclo = pack(iclo, iclo /= -1)
!             rclo = pack(rclo, rclo /= -1)
!             dclo = pack(dclo, dclo >= 0.0d0)
!
!             !---------------------------------------------------------------------
!             ! Update isma, iclo and rclo with indices of remote parcels:
!             do m = 1, n_mix
!
!                 m = m + 1
!
!                 is = isma(n)
!                 ic = iclo(n)
!
!                 if ((is > source%local_num) .and. (ic > source%local_num)) then
!                     ! A remote small parcel points to another remote small parcel. The remotes do not necessarily
!                     ! need to be the same. Also, it can be a dual-link. As the same distance is evaluated on
!                     ! the other two MPI ranks, *this* MPI rank must set the distance between the parcels to the
!                     ! maximum value as otherwise the function "find_closest_parcel_globally" may think the
!                     ! parcels belong to *this* MPI rank due to round-offs in the distance calculation. The
!                     ! function "find_closest_parcel_globally" always sets "rclo" to the MPI source.
!                     dclo(m) = huge(0.0d0) ! huge(x) returns the maximum value of this type
!                 else if (ic > source%local_num) then
!                     ! A local small parcel points to a remote small parcel.
!                     ! The index *ic* is larger than the local number of parcels
!                     ! we must update *iclo(m)* and *rclo(m)* with the index of the parcel stored
!                     ! on the other MPI rank.
!                     iclo(m) = pidsma(ic)
!                     rclo(m) = rsma(ic)
!                 endif
!             enddo
!
!             !------------------------------------------------------------------
!             ! Sanity check: Indices of small parcels must be smaller equal to the
!             ! number of local parcels:
!             do m = 1, n_local_mix
!                 if (isma(m) > source%local_num) then
!                     call mpi_exit_on_error(&
!                         'in in parcel_mixing::find_locally: Small parcel index out of range.')
!                 endif
!             enddo

        end subroutine find_locally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine find_globally
            print *, "TODO"
        end subroutine find_globally

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine actual_mixing(source, dest, nmix)
            class(pc_type),  intent(inout) :: source
            class(pc_type),  intent(inout) :: dest
            integer,         intent(in)    :: nmix
            integer                        :: l, m, n, ic, is
            double precision               :: buoym(nmix), vortm(3, nmix), vm(nmix), vmix
            integer                        :: lclo(dest%local_num)

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
