! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_ellipse_merge
    use parcel_nearest
    use constants, only : pi, zero, one, two, four
    use parcel_ellipse, only : ellipse_pc_type
    use options, only : parcel, verbose
    use parcel_bc
    use parcel_ops, only : get_delx_across_periodic   &
                         , get_dely_across_periodic
!     use timer, only : start_timer, stop_timer
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use parcel_mpi, only : parcel_communicate
    implicit none

    private

    integer, allocatable :: loca(:)

    public :: surface_parcel_merge

    contains

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        subroutine surface_parcel_merge(spar)
            type(ellipse_pc_type), intent(inout) :: spar
            integer, allocatable, dimension(:)   :: isma
            integer, allocatable, dimension(:)   :: iclo
            integer, allocatable, dimension(:)   :: inva
            integer                              :: n_merge ! number of merges
            integer                              :: n_invalid
#if defined (ENABLE_VERBOSE) && !defined (NDEBUG)
            integer                              :: orig_num

            orig_num = spar%total_num
#endif

            ! find parcels to merge
            call find_nearest(spar, isma, iclo, inva, n_merge, n_invalid)

!             call start_timer(merge_timer)

!             n_parcel_merges = n_parcel_merges + n_merge

            if (n_merge > 0) then
                allocate(loca(spar%local_num))
!                 call collect_merge_stats(iclo, n_merge)
            endif


            if (n_merge > 0) then
                ! merge small parcels into other parcels
                call group_merge(spar, isma, iclo, n_merge)
            endif

            if (n_merge + n_invalid > 0) then
                ! overwrite all small and invalid parcels -- all small parcels are now invalid too
                call spar%delete(inva, n_merge + n_invalid)
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
            spar%total_num = spar%local_num
            call mpi_blocking_reduce(spar%total_num, MPI_SUM, world)

#if defined (ENABLE_VERBOSE) && !defined (NDEBUG)
            if (verbose .and. (world%rank == world%root)) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", orig_num,  &
                      "...", spar%total_num
            endif
#endif
            call parcel_communicate(spar)

!             call stop_timer(merge_timer)

        end subroutine surface_parcel_merge

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Geometric merge.
        ! @param[inout] spar is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine group_merge(spar, isma, iclo, n_merge)
            type(ellipse_pc_type), intent(inout) :: spar
            integer,               intent(in)    :: isma(0:)
            integer,               intent(in)    :: iclo(:)
            integer,               intent(in)    :: n_merge
            integer                              :: m, ic, is, l, n
            double precision                     :: factor
            double precision                     :: x0(n_merge), y0(n_merge)
            double precision                     :: posm(2, n_merge), delx, amerge, dely, mu
            double precision                     :: buoym(n_merge), vortm(2, n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                     :: hum(n_merge)
#endif
            double precision                     :: Bm(3, n_merge), am(n_merge), vm(n_merge), cm(n_merge)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                ! This is different to the serial version: We must apply a periodic
                ! shift as a small parcel might be across a periodic boundary.
                is = isma(m)
                call apply_periodic_bc(spar%position(:, is))

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! am will contain the total area of the merged parcel
                    am(l) = spar%area(ic)

                    ! cm will contain the total circulation of the merged parcel
                    cm(l) = spar%circ(ic)

                    vm(l) = spar%volume(ic)


                    !x0 stores the x centre of the other parcel
                    x0(l) = spar%position(1, ic)

                    !y0 stores the y centre of the other parcel
                    y0(l) = spar%position(2, ic)

                    ! posm(:, l) will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    posm(:, l) = zero

                    ! buoyancy and humidity
                    buoym(l) = spar%area(ic) * spar%buoyancy(ic)
#ifndef ENABLE_DRY_MODE
                    hum(l) = spar%area(ic) * spar%humidity(ic)
#endif
                    vortm(:, l) = spar%area(ic) * spar%vorticity(:, ic)

                    Bm(:, l) = zero
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                am(n) = am(n) + spar%area(is) !Accumulate area of merged parcel
                cm(n) = cm(n) + spar%circ(is)

                vm(n) = vm(n) + spar%volume(is)

                ! works across periodic edge
                delx = get_delx_across_periodic(spar%position(1, is), x0(n))
                dely = get_dely_across_periodic(spar%position(2, is), y0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic)) and v(is)*(y(is)-y(ic))
                posm(1, n) = posm(1, n) + spar%area(is) * delx
                posm(2, n) = posm(2, n) + spar%area(is) * dely

                ! Accumulate buoyancy and humidity
                buoym(n) = buoym(n) + spar%area(is) * spar%buoyancy(is)
#ifndef ENABLE_DRY_MODE
                hum(n) = hum(n) + spar%area(is) * spar%humidity(is)
#endif
                vortm(:, n) = vortm(:, n) + spar%area(is) * spar%vorticity(:, is)
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / am(m)
                amerge = one / am(m)

                ! need to sanitise input and output, but first to determine input
                posm(:, m) = - amerge * posm(:, m)

                call apply_periodic_bc(posm(:, m))

                ! x and y centre of merged parcel, modulo periodicity
                posm(1, m) = get_delx_across_periodic(x0(m), posm(1, m))
                posm(2, m) = get_dely_across_periodic(y0(m), posm(2, m))

                ! need to correct position
                call apply_periodic_bc(posm(:, m))

                ! buoyancy and humidity
                buoym(m) = amerge * buoym(m)
#ifndef ENABLE_DRY_MODE
                hum(m) = amerge * hum(m)
#endif
                vortm(:, m) = amerge * vortm(:, m)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    amerge = one / am(l)

                    delx = get_delx_across_periodic(spar%position(1, ic), posm(1, l))
                    dely = get_dely_across_periodic(spar%position(2, ic), posm(2, l))

                    mu = spar%area(ic) * amerge
                    Bm(1, l) = mu * (four * delx ** 2   + spar%B(1, ic))
                    Bm(2, l) = mu * (four * delx * dely + spar%B(2, ic))
                    Bm(3, l) = mu * (four * dely ** 2   + spar%B(3, ic))

                    spar%area(ic) = am(l)
                    spar%circ(ic) = cm(l)
                    spar%volume(ic) = vm(l)
                    spar%position(1, ic) = posm(1, l)
                    spar%position(2, ic) = posm(2, l)

                    spar%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    spar%humidity(ic) = hum(l)
#endif
                    spar%vorticity(:, ic) = vortm(:, l)

                endif

                is = isma(m)
                n = loca(ic)

                amerge = one / am(n)

                delx = get_delx_across_periodic(spar%position(1, is), posm(1, n))
                dely = get_dely_across_periodic(spar%position(2, is), posm(2, n))

                ! area fraction A_{is} / A
                mu = amerge * spar%area(is)

                Bm(1, n) = Bm(1, n) + mu * (four * delx ** 2   + spar%B(1, is))
                Bm(2, n) = Bm(2, n) + mu * (four * delx * dely + spar%B(2, is))
                Bm(3, n) = Bm(3, n) + mu * (four * dely ** 2   + spar%B(3, is))
            enddo

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! normalize such that determinant of the merger is (ab)**2
                    ! ab / sqrt(det(B))
                    factor = spar%get_ab(am(l)) / dsqrt(dabs(Bm(1, l) * Bm(3, l) - Bm(2, l) ** 2))

                    spar%B(:, ic) = Bm(:, l) * factor
                endif
            enddo

        end subroutine group_merge

end module parcel_ellipse_merge
