! =============================================================================
!                       Module to merge ellipsoids
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_merging
    use parcel_nearest
    use parameters, only : vmin
    use constants, only : pi, zero, one, two, five, f13
    use parcel_container, only : parcels                    &
                               , n_parcels                  &
                               , n_total_parcels            &
                               , parcel_replace             &
                               , get_delx_across_periodic   &
                               , get_dely_across_periodic   &
                               , parcel_delete
    use parcel_ellipsoid, only : get_abc
    use options, only : parcel
    use datatypes, only : int64
#if defined (ENABLE_VERBOSE) && !defined (NDEBUG)
    use options, only : verbose
#endif
    use parcel_bc, only : apply_periodic_bc
    use parcel_mpi, only : parcel_communicate
    use mpi_timer, only : start_timer, stop_timer
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    integer :: merge_timer

    ! number of parcel merges (is reset in every write step)
    integer(kind=int64) :: n_parcel_merges = 0

    ! number of merging parcels (up to 7 supported, all others are put into index 7)
    ! note that array index 1 corresponds to 2-way merging
    integer(kind=int64) :: n_way_parcel_mergers(7) = 0

    ! number of big iclo neighbours (number of small is n_merge - n_big_close)
    integer :: n_big_close = 0

    integer, allocatable :: loca(:)

    private :: geometric_merge,     &
               do_group_merge,      &
               collect_merge_stats, &
               loca

    contains

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        subroutine parcel_merge
            integer, allocatable, dimension(:) :: isma
            integer, allocatable, dimension(:) :: iclo
            integer, allocatable, dimension(:) :: inva
            integer                            :: n_merge ! number of merges
            integer                            :: n_invalid
#if defined (ENABLE_VERBOSE) && !defined (NDEBUG)
            integer(kind=int64)                :: orig_num

            orig_num = n_total_parcels
#endif

            ! find parcels to merge
            call find_nearest(isma, iclo, inva, n_merge, n_invalid)

            call start_timer(merge_timer)

            n_parcel_merges = n_parcel_merges + n_merge

            if (n_merge > 0) then
                allocate(loca(n_parcels))
                call collect_merge_stats(iclo, n_merge)
            endif


            if (n_merge > 0) then
                ! merge small parcels into other parcels
                call geometric_merge(isma, iclo, n_merge)
            endif

            if (n_merge + n_invalid > 0) then
                ! overwrite all small and invalid parcels -- all small parcels are now invalid too
                call parcel_delete(inva, n_merge + n_invalid)
                deallocate(inva)
            endif

            if (allocated(isma)) then
                deallocate(isma)
                deallocate(iclo)
            endif

            if (allocated(loca)) then
                deallocate(loca)
            endif

            ! After this operation the root MPI process knows the new
            ! number of parcels in the simulation
            n_total_parcels = n_parcels
            call mpi_blocking_reduce(n_total_parcels, MPI_SUM, world)

#if defined (ENABLE_VERBOSE) && !defined (NDEBUG)
            if (verbose .and. (world%rank == world%root)) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", orig_num,  &
                      "...", n_total_parcels
            endif
#endif
            call parcel_communicate

            call stop_timer(merge_timer)

        end subroutine parcel_merge

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Actual merge.
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        ! @param[out] Bm are the B matrix entries of the mergers
        ! @param[out] vm are the volumes of the mergers
        subroutine do_group_merge(isma, iclo, n_merge, Bm, vm)
            integer,          intent(in)    :: isma(0:)
            integer,          intent(in)    :: iclo(:)
            integer,          intent(in)    :: n_merge
            integer                         :: m, ic, is, l, n
            double precision                :: x0(n_merge), y0(n_merge)
            double precision                :: posm(3, n_merge)
            double precision                :: delx, vmerge, dely, delz, mu
            double precision                :: buoym(n_merge), vortm(3, n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                :: hum(n_merge)
#endif
#ifdef ENABLE_LABELS
            double precision                :: labelm(n_merge)
            double precision                :: dilm(n_merge)
            double precision                :: rn
#endif
            double precision, intent(out)   :: Bm(6, n_merge) ! B11, B12, B13, B22, B23, B33
            double precision, intent(out)   :: vm(n_merge)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                ! This is different to the serial version: We must apply a periodic
                ! shift as a small parcel might be across a periodic boundary.
                is = isma(m)
                call apply_periodic_bc(parcels%position(:, is))

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! vm will contain the total volume of the merged parcel
                    vm(l) = parcels%volume(ic)

                    !x0 stores the x centre of the other parcel
                    x0(l) = parcels%position(1, ic)

                    !y0 stores the y centre of the other parcel
                    y0(l) = parcels%position(2, ic)

                    ! posm(1, :) will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    posm(1, l) = zero

                    ! posm(2, :) will sum v(is)*(y(is)-y(ic)) modulo periodicity
                    posm(2, l) = zero

                    ! posm(3, :) will contain v(ic)*z(ic)+sum{v(is)*z(is)}
                    posm(3, l) = parcels%volume(ic) * parcels%position(3, ic)

                    ! buoyancy and humidity
                    buoym(l) = parcels%volume(ic) * parcels%buoyancy(ic)
#ifndef ENABLE_DRY_MODE
                    hum(l) = parcels%volume(ic) * parcels%humidity(ic)
#endif
                    vortm(:, l) = parcels%volume(ic) * parcels%vorticity(:, ic)

                    Bm(:, l) = zero

#ifdef ENABLE_LABELS
                    labelm(l) = parcels%label(ic)

                    dilm(l) = parcels%dilution(ic)
#endif
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel

#ifdef ENABLE_LABELS
                ! Dilute the parcel when volume is added for now
                ! This could be optimised moving it to later in code
                call random_number(rn) 
                if(rn*(vm(n)+parcels%volume(is))<vm(n)) then
                   dilm(n)=dilm(n)+log(vm(n)/(vm(n)+parcels%volume(is)))
                else
                   labelm(n)=parcels%label(is)
                   dilm(n)=parcels%dilution(is)+log(parcels%volume(is)/(vm(n)+parcels%volume(is)))
                endif
#endif
                vm(n) = vm(n) + parcels%volume(is) !Accumulate volume of merged parcel

                ! works across periodic edge
                delx = get_delx_across_periodic(parcels%position(1, is), x0(n))
                dely = get_dely_across_periodic(parcels%position(2, is), y0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic)) and v(is)*(y(is)-y(ic))
                posm(1, n) = posm(1, n) + parcels%volume(is) * delx
                posm(2, n) = posm(2, n) + parcels%volume(is) * dely

                ! Accumulate v(ic)*z(ic)+sum{v(is)*z(is)}
                posm(3, n) = posm(3, n) + parcels%volume(is) * parcels%position(3, is)

                ! Accumulate buoyancy and humidity
                buoym(n) = buoym(n) + parcels%volume(is) * parcels%buoyancy(is)
#ifndef ENABLE_DRY_MODE
                hum(n) = hum(n) + parcels%volume(is) * parcels%humidity(is)
#endif
                vortm(:, n) = vortm(:, n) + parcels%volume(is) * parcels%vorticity(:, is)
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmerge = one / vm(m)

                ! need to sanitise input and output, but first to determine input
                posm(1, m) = - vmerge * posm(1, m)
                posm(2, m) = - vmerge * posm(2, m)

                call apply_periodic_bc(posm(:, m))
                ! x and y centre of merged parcel, modulo periodicity
                posm(1, m) = get_delx_across_periodic(x0(m), posm(1, m))
                posm(2, m) = get_dely_across_periodic(y0(m), posm(2, m))

                ! z centre of merged parcel
                posm(3, m) = vmerge * posm(3, m)

                ! need to correct position
                call apply_periodic_bc(posm(:, m))

                ! buoyancy and humidity
                buoym(m) = vmerge * buoym(m)
#ifndef ENABLE_DRY_MODE
                hum(m) = vmerge * hum(m)
#endif
                vortm(:, m) = vmerge * vortm(:, m)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    vmerge = one / vm(l)

                    delx = get_delx_across_periodic(parcels%position(1, ic), posm(1, l))
                    dely = get_dely_across_periodic(parcels%position(2, ic), posm(2, l))
                    delz = parcels%position(3, ic) - posm(3, l)

                    mu = parcels%volume(ic) * vmerge


                    Bm(1, l) = mu * (five * delx ** 2   + parcels%B(1, ic))
                    Bm(2, l) = mu * (five * delx * dely + parcels%B(2, ic))
                    Bm(3, l) = mu * (five * delx * delz + parcels%B(3, ic))
                    Bm(4, l) = mu * (five * dely ** 2   + parcels%B(4, ic))
                    Bm(5, l) = mu * (five * dely * delz + parcels%B(5, ic))
                    Bm(6, l) = mu * (five * delz ** 2   + parcels%B(6, ic))

                    parcels%volume(ic)  = vm(l)
                    parcels%position(1, ic) = posm(1, l)
                    parcels%position(2, ic) = posm(2, l)
                    parcels%position(3, ic) = posm(3, l)

                    parcels%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ic) = hum(l)
#endif
#ifdef ENABLE_LABELS
                    parcels%label(ic) = labelm(l)
                    parcels%dilution(ic) = dilm(l)
#endif
                    parcels%vorticity(:, ic) = vortm(:, l)

                endif

                is = isma(m)
                n = loca(ic)

                vmerge = one / vm(n)

                delx = get_delx_across_periodic(parcels%position(1, is), posm(1, n))
                dely = get_dely_across_periodic(parcels%position(2, is), posm(2, n))
                delz = parcels%position(3, is) - posm(3, n)

                ! volume fraction V_{is} / V
                mu = vmerge * parcels%volume(is)

                Bm(1, n) = Bm(1, n) + mu * (five * delx ** 2   + parcels%B(1, is))
                Bm(2, n) = Bm(2, n) + mu * (five * delx * dely + parcels%B(2, is))
                Bm(3, n) = Bm(3, n) + mu * (five * delx * delz + parcels%B(3, is))
                Bm(4, n) = Bm(4, n) + mu * (five * dely ** 2   + parcels%B(4, is))
                Bm(5, n) = Bm(5, n) + mu * (five * dely * delz + parcels%B(5, is))
                Bm(6, n) = Bm(6, n) + mu * (five * delz ** 2   + parcels%B(6, is))
            enddo

        end subroutine do_group_merge

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Geometric merging -- called by subroutine merge_parcels.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine geometric_merge(isma, iclo, n_merge)
            integer,         intent(in) :: isma(0:)
            integer,         intent(in) :: iclo(:)
            integer,         intent(in) :: n_merge
            integer                     :: m, ic, l
            double precision            :: factor, detB
            double precision            :: B(6, n_merge), &
                                        V(n_merge)

            call do_group_merge(isma, iclo, n_merge, B, V)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! normalize such that determinant of the merger is (abc)**2
                    ! ((abc)**2 / det(B))^(1/3)
                    detB = B(1, l) * (B(4, l) * B(6, l) - B(5, l) ** 2) &
                         - B(2, l) * (B(2, l) * B(6, l) - B(3, l) * B(5, l)) &
                         + B(3, l) * (B(2, l) * B(5, l) - B(3, l) * B(4, l))

                    factor = (get_abc(V(l)) ** 2 / detB) ** f13

                    parcels%B(:, ic) = B(:, l) * factor
                endif
            enddo

        end subroutine geometric_merge

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine collect_merge_stats(iclo, n_merge)
            integer, allocatable, dimension(:) :: iclo
            integer                            :: n_merge
            integer                            :: ic, m, j, n_count

            loca = 0

            !------------------------------------------------------------------
            ! Find unique 'iclo' indices and the total number of parcels
            ! merging with them:
            do m = 1, n_merge
                ic = iclo(m)
                n_big_close = n_big_close + merge(1, 0, parcels%volume(ic) > vmin)
                loca(ic) = loca(ic) + 1
            enddo

            !------------------------------------------------------------------
            ! Count the number of 2-, 3-, 4- etc way merging:
            do m = 1, n_merge
                ic = iclo(m)
                n_count = loca(ic)
                ! all mergers involving more than size(n_way_parcel_mergers) parcels are added together
                if (n_count > 0) then
                    loca(ic) = -1
                    j = min(size(n_way_parcel_mergers), n_count)
                    n_way_parcel_mergers(j) = n_way_parcel_mergers(j) + 1
                endif
            enddo

        end subroutine collect_merge_stats

end module parcel_merging
