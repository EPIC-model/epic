! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_merge
    use parcel_nearest
    use constants, only : pi, zero, one, two, four
    use parcel_container, only : get_delx
    use dynamic_parcels, only : n_parcels, parcels
    use parcel_ellipsoid, only : ellipsoid_parcel_type
    use parcel_ellipse, only : get_B22, get_ab
    use options, only : parcel, verbose
    use parcel_bc
    use timer, only : start_timer, stop_timer
    implicit none

    integer:: merge_timer

    private :: geometric_merge, &
               do_group_merge,  &
               pack_parcels

    contains

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        ! @param[inout] parcels is the parcel container
        subroutine merge_ellipses(parcels)
            class(ellipsoid_parcel_type), intent(inout) :: parcels
            integer, allocatable, dimension(:)         :: isma
            integer, allocatable, dimension(:)         :: iclo
            integer                                    :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(isma, iclo, n_merge)

            call start_timer(merge_timer)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", n_parcels, &
                      "...", n_parcels - n_merge
            endif
#endif

            if (n_merge > 0) then
                ! merge small parcels into other parcels
                call geometric_merge(parcels, isma, iclo, n_merge)

                ! overwrite invalid parcels
                call pack_parcels(isma, n_merge)
            endif

            if (allocated(isma)) then
                deallocate(isma)
                deallocate(iclo)
            endif

            call stop_timer(merge_timer)

        end subroutine merge_ellipses


        ! Actual merge.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine do_group_merge(parcels, isma, iclo, n_merge)
            class(ellipsoid_parcel_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, is, l, n
            integer                                    :: loca(n_parcels)
            double precision                           :: x0(n_merge)
            double precision                           :: delx, vmerge, dely, B22, mu

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! vm will contain the total volume of the merged parcel
                    parcels%merge_volume(l) = parcels%volume(ic)

                    !x0 stores the x centre of the other parcel
                    x0(l) = parcels%position(1, ic)

                    ! posm(1, l) will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    parcels%merge_position(1, l) = zero

                    ! posm(2, l) will contain v(ic)*z(ic)+sum{v(is)*z(is)}
                    parcels%merge_position(2, l) = parcels%volume(ic) * parcels%position(2, ic)

                    call parcels%assign_p2m(ic, l, parcels%volume(ic))

                    parcels%merge_B(1, l) = zero
                    parcels%merge_B(2, l) = zero
                    parcels%merge_B(3, l) = zero !B22
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                parcels%merge_volume(n) = parcels%merge_volume(n) + parcels%volume(is) !Accumulate volume of merged parcel

                ! works across periodic edge
                delx = get_delx(parcels%position(1, is), x0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic))
                parcels%merge_position(1, n) = parcels%merge_position(1, n) + parcels%volume(is) * delx

                ! Accumulate v(ic)*z(ic)+sum{v(is)*z(is)}
                parcels%merge_position(2, n) = parcels%merge_position(2, n) + parcels%volume(is) * parcels%position(2, is)

                call parcels%add_p2m(is, n, parcels%volume(is))
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmerge = one / parcels%merge_volume(m)

                ! x centre of merged parcel, modulo periodicity
                parcels%merge_position(1, m) = vmerge * parcels%merge_position(1, m)

                parcels%merge_position(1, m) = x0(m) + parcels%merge_position(1, m)

                ! z centre of merged parcel
                parcels%merge_position(2, m) = vmerge * parcels%merge_position(2, m)

                ! need to correct position
                call apply_periodic_bc(parcels%merge_position(:, m))

                ! normalise other properties by vmerge
                call parcels%assign_m2m(m, vmerge)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    vmerge = one / parcels%merge_volume(l)

                    B22 = get_B22(parcels%B(1, ic), parcels%B(2, ic), parcels%volume(ic))

                    delx = get_delx(parcels%position(1, ic), parcels%merge_position(1, l))
                    dely = parcels%position(2, ic) - parcels%merge_position(2, l)

                    mu = parcels%volume(ic) * vmerge
                    parcels%merge_B(1, l) = mu * (four * delx ** 2 + parcels%B(1, ic))
                    parcels%merge_B(2, l) = mu * (four * delx * dely + parcels%B(2, ic))
                    parcels%merge_B(3, l) = mu * (four * dely ** 2 + B22)

                    parcels%volume(ic)  = parcels%merge_volume(l)
                    parcels%position(1, ic) = parcels%merge_position(1, l)
                    parcels%position(2, ic) = parcels%merge_position(2, l)

                    call parcels%assign_m2p(l, ic)
                endif

                is = isma(m)
                n = loca(ic)

                vmerge = one / parcels%merge_volume(n)

                delx = get_delx(parcels%position(1, is), parcels%merge_position(1, n))
                dely = parcels%position(2, is) - parcels%merge_position(2, n)

                B22 = get_B22(parcels%B(1, is), parcels%B(2, is), parcels%volume(is))

                ! volume fraction A_{is} / A
                mu = vmerge * parcels%volume(is)

                parcels%merge_B(1, n) = parcels%merge_B(1, n) + mu * (four * delx ** 2   + parcels%B(1, is))
                parcels%merge_B(2, n) = parcels%merge_B(2, n) + mu * (four * delx * dely + parcels%B(2, is))
                parcels%merge_B(3, n) = parcels%merge_B(3, n) + mu * (four * dely ** 2   + B22)
            enddo

        end subroutine do_group_merge


        ! Geometric merging -- called by subroutine merge_ellipses.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine geometric_merge(parcels, isma, iclo, n_merge)
            class(ellipsoid_parcel_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, l
            integer                                    :: loca(n_parcels)
            double precision                           :: factor

            call parcels%merge_alloc(n_merge)

            call do_group_merge(parcels, isma, iclo, n_merge)

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
                    factor = get_ab(parcels%merge_volume(l)) / &
                             sqrt(parcels%merge_B(1, l) * parcels%merge_B(3, l) - parcels%merge_B(2, l) ** 2)

                    parcels%B(1, ic) = parcels%merge_B(1, l) * factor
                    parcels%B(2, ic) = parcels%merge_B(2, l) * factor

                    call apply_periodic_bc(parcels%position(:, ic))
                endif
            enddo

            call parcels%merge_dealloc

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
        subroutine pack_parcels(isma, n_merge)
            integer, intent(in) :: isma(0:)
            integer, intent(in) :: n_merge
            integer             :: k, l, m

            ! l points always to the last valid parcel
            l = n_parcels

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
                call parcels%replace(isma(m), l)

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
            n_parcels = n_parcels - n_merge

        end subroutine pack_parcels


end module parcel_merge
