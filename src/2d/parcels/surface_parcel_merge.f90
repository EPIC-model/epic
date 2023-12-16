! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
! =============================================================================
module surface_parcel_merge_mod
    use surface_parcel_nearest
    use constants, only : pi, zero, one, two, four
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_top_parcels, n_bot_parcels   &
                                       , top_parcels, bot_parcels       &
                                       , surface_parcel_replace
    use parcel_ops, only : get_delx
    use options, only : parcel, verbose
    use surface_parcel_bc
    implicit none

    private :: geometric_merge, &
               do_group_merge,  &
               pack_parcels

    contains

        subroutine surface_parcel_merge

            call merge_lines(n_bot_parcels, bot_parcels)
            call merge_lines(n_top_parcels, top_parcels)

        end subroutine surface_parcel_merge

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        ! @param[inout] parcels is the parcel container
        subroutine merge_lines(n_par, spar)
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer, allocatable, dimension(:)                 :: isma
            integer, allocatable, dimension(:)                 :: iclo
            integer                                            :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(n_par, spar, isma, iclo, n_merge)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)",                                   &
                      "no. surface parcels before and after merge: ", n_par, &
                      "...", n_par - n_merge
            endif
#endif

            if (n_merge > 0) then
                ! merge small parcels into other parcels
                call geometric_merge(n_par, spar, isma, iclo, n_merge)

                ! overwrite invalid parcels
                call pack_parcels(n_par, spar, isma, n_merge)
            endif

            if (allocated(isma)) then
                deallocate(isma)
                deallocate(iclo)
            endif

        end subroutine merge_lines


        ! Actual merge.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        ! @param[out] lm are the lengths of the mergers
        subroutine do_group_merge(n_par, spar, isma, iclo, n_merge, lm)
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(0:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: n_merge
            double precision,                    intent(out)   :: lm(n_merge)
            integer                                            :: m, ic, is, l, n
            integer                                            :: loca(n_par)
            double precision                                   :: x0(n_merge)
            double precision                                   :: posm(n_merge), delx, lmerge
            double precision                                   :: buoym(n_merge), vortm(n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                                   :: hum(n_merge)
#endif

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! lm will contain the total length of the merged parcel
                    lm(l) = spar%length(ic)

                    !x0 stores the x centre of the other parcel
                    x0(l) = spar%position(ic)

                    ! posm(l) will sum length(is)*(x(is)-x(ic)) modulo periodicity
                    posm(l) = zero

                    ! buoyancy and humidity
                    buoym(l) = spar%length(ic) * spar%buoyancy(ic)
#ifndef ENABLE_DRY_MODE
                    hum(l) = spar%length(ic) * spar%humidity(ic)
#endif
                    vortm(l) = spar%length(ic) * spar%vorticity(ic)

                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                lm(n) = lm(n) + spar%length(is) !Accumulate length of merged parcel

                ! works across periodic edge
                delx = get_delx(spar%position(is), x0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic))
                posm( n) = posm(n) + spar%length(is) * delx

                ! Accumulate buoyancy and humidity
                buoym(n) = buoym(n) + spar%length(is) * spar%buoyancy(is)
#ifndef ENABLE_DRY_MODE
                hum(n) = hum(n) + spar%length(is) * spar%humidity(is)
#endif
                vortm(n) = vortm(n) + spar%length(is) * spar%vorticity(is)
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / lm(m)
                lmerge = one / lm(m)

                ! need to sanitise input and output, but first to determine input
                posm(m) = - lmerge * posm(m)

                call apply_surface_periodic_bc(posm(m))

                ! x centre of merged parcel, modulo periodicity
                posm(m) = get_delx(x0(m), posm(m))

                ! need to correct position
                call apply_surface_periodic_bc(posm(m))

                ! buoyancy and humidity
                buoym(m) = lmerge * buoym(m)
#ifndef ENABLE_DRY_MODE
                hum(m) = lmerge * hum(m)
#endif
                vortm(m) = lmerge * vortm(m)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    delx = get_delx(spar%position(ic), posm(l))

                    spar%length(ic)  = lm(l)
                    spar%position(ic) = posm(l)

                    spar%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    spar%humidity(ic) = hum(l)
#endif
                    spar%vorticity(ic) = vortm(l)

                endif

                is = isma(m)
                n = loca(ic)

                delx = get_delx(spar%position(is), posm(n))

            enddo

        end subroutine do_group_merge


        ! Geometric merging -- called by subroutine merge_ellipses.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine geometric_merge(n_par, spar, isma, iclo, n_merge)
            integer,                             intent(in)    :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(0:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: n_merge
            integer                                            :: m, ic, l
            integer                                            :: loca(n_par)
            double precision                                   :: length(n_merge)

            call do_group_merge(n_par, spar, isma, iclo, n_merge, length)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    call apply_surface_periodic_bc(spar%position(ic))
                endif
            enddo

        end subroutine geometric_merge


        ! This algorithm replaces invalid parcels with valid parcels
        ! from the end of the container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] n_merge is the array size of isma and iclo
        ! @pre
        !   - isma must be sorted in ascending order
        !   - isma must be contiguously filled
        !   The above preconditions must be fulfilled so that the
        !   parcel pack algorithm works correctly.
        subroutine pack_parcels(n_par, spar, isma, n_merge)
            integer,                             intent(inout) :: n_par
            type(surface_parcel_container_type), intent(inout) :: spar
            integer,                             intent(in)    :: isma(0:)
            integer,                             intent(in)    :: n_merge
            integer                                            :: k, l, m

            ! l points always to the last valid parcel
            l = n_par

            ! k points always to last invalid parcel in isma
            k = n_merge

            ! find last parcel which is not invalid
            do while ((k > 0) .and. (l == isma(k)))
                l = l - 1
                k = k - 1
            enddo

            if (l == 0) then
                print *, "Error: All parcels are invalid."
                stop
            endif

            ! replace invalid parcels with the last valid parcel
            m = 1

            do while (m <= k)
                ! invalid parcel; overwrite *isma(m)* with last valid parcel *l*
                call surface_parcel_replace(isma(m), l, spar)

                l = l - 1

                ! find next valid last parcel
                do while ((k > 0) .and. (l == isma(k)))
                    l = l - 1
                    k = k - 1
                enddo

                ! next invalid
                m = m + 1
            enddo

            ! update number of valid parcels
            n_par = n_par - n_merge

        end subroutine pack_parcels


end module surface_parcel_merge_mod
