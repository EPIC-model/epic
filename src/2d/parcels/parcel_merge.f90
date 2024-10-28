! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_merge
    use parcel_nearest
    use constants, only : pi, zero, one, two, four
    use parcel_container, only : parcel_container_type  &
                               , n_parcels              &
                               , parcel_replace         &
                               , get_delx
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
            type(parcel_container_type), intent(inout) :: parcels
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
        ! @param[out] B11m are the B11 matrix entries of the mergers
        ! @param[out] B12m are the B12 matrix entries of the mergers
        ! @param[out] B22m are the B22 matrix entries of the mergers
        ! @param[out] vm are the volumes of the mergers
        subroutine do_group_merge(parcels, isma, iclo, n_merge, B11m, B12m, B22m, vm)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, is, l, n
            integer                                    :: loca(n_parcels)
            double precision                           :: x0(n_merge)
            double precision                           :: posm(2, n_merge), delx, vmerge, dely, B22, mu
            double precision                           :: buoym(n_merge), vortm(n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                           :: hum(n_merge)
#endif
            double precision,            intent(out)   :: B11m(n_merge), B12m(n_merge), B22m(n_merge), &
                                                          vm(n_merge)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! vm will contain the total volume of the merged parcel
                    vm(l) = parcels%volume(ic)

                    !x0 stores the x centre of the other parcel
                    x0(l) = parcels%position(1, ic)

                    ! posm(1, l) will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    posm(1, l) = zero

                    ! posm(2, l) will contain v(ic)*z(ic)+sum{v(is)*z(is)}
                    posm(2, l) = parcels%volume(ic) * parcels%position(2, ic)

                    ! buoyancy and humidity
                    buoym(l) = parcels%volume(ic) * parcels%buoyancy(ic)
#ifndef ENABLE_DRY_MODE
                    hum(l) = parcels%volume(ic) * parcels%humidity(ic)
#endif
                    vortm(l) = parcels%volume(ic) * parcels%vorticity(ic)

                    B11m(l) = zero
                    B12m(l) = zero
                    B22m(l) = zero
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                vm(n) = vm(n) + parcels%volume(is) !Accumulate volume of merged parcel

                ! works across periodic edge
                delx = get_delx(parcels%position(1, is), x0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic))
                posm(1, n) = posm(1, n) + parcels%volume(is) * delx

                ! Accumulate v(ic)*z(ic)+sum{v(is)*z(is)}
                posm(2, n) = posm(2, n) + parcels%volume(is) * parcels%position(2, is)

                ! Accumulate buoyancy and humidity
                buoym(n) = buoym(n) + parcels%volume(is) * parcels%buoyancy(is)
#ifndef ENABLE_DRY_MODE
                hum(n) = hum(n) + parcels%volume(is) * parcels%humidity(is)
#endif
                vortm(n) = vortm(n) + parcels%volume(is) * parcels%vorticity(is)
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmerge = one / vm(m)

                ! x centre of merged parcel, modulo periodicity
                posm(1, m) = vmerge * posm(1, m)

                posm(1, m) = x0(m) + posm(1, m)

                ! z centre of merged parcel
                posm(2, m) = vmerge * posm(2, m)

                ! need to correct position
                call apply_periodic_bc(posm(:, m))

                ! buoyancy and humidity
                buoym(m) = vmerge * buoym(m)
#ifndef ENABLE_DRY_MODE
                hum(m) = vmerge * hum(m)
#endif
                vortm(m) = vmerge * vortm(m)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    vmerge = one / vm(l)

                    B22 = get_B22(parcels%B(1, ic), parcels%B(2, ic), parcels%volume(ic))

                    delx = get_delx(parcels%position(1, ic), posm(1, l))
                    dely = parcels%position(2, ic) - posm(2, l)

                    mu = parcels%volume(ic) * vmerge
                    B11m(l) = mu * (four * delx ** 2 + parcels%B(1, ic))
                    B12m(l) = mu * (four * delx * dely + parcels%B(2, ic))
                    B22m(l) = mu * (four * dely ** 2 + B22)

                    parcels%volume(ic)  = vm(l)
                    parcels%position(1, ic) = posm(1, l)
                    parcels%position(2, ic) = posm(2, l)

                    parcels%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ic) = hum(l)
#endif
                    parcels%vorticity(ic) = vortm(l)

                endif

                is = isma(m)
                n = loca(ic)

                vmerge = one / vm(n)

                delx = get_delx(parcels%position(1, is), posm(1, n))
                dely = parcels%position(2, is) - posm(2, n)

                B22 = get_B22(parcels%B(1, is), parcels%B(2, is), parcels%volume(is))

                ! volume fraction A_{is} / A
                mu = vmerge * parcels%volume(is)

                B11m(n) = B11m(n) + mu * (four * delx ** 2   + parcels%B(1, is))
                B12m(n) = B12m(n) + mu * (four * delx * dely + parcels%B(2, is))
                B22m(n) = B22m(n) + mu * (four * dely ** 2   + B22)
            enddo

        end subroutine do_group_merge


        ! Geometric merging -- called by subroutine merge_ellipses.
        ! @param[inout] parcels is the parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine geometric_merge(parcels, isma, iclo, n_merge)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, l
            integer                                    :: loca(n_parcels)
            double precision                           :: factor
            double precision                           :: B11(n_merge), &
                                                          B12(n_merge), &
                                                          B22(n_merge), &
                                                          V(n_merge)

            call do_group_merge(parcels, isma, iclo, n_merge, B11, B12, B22, V)

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
                    factor = get_ab(V(l)) / sqrt(B11(l) * B22(l) - B12(l) ** 2)

                    parcels%B(1, ic) = B11(l) * factor
                    parcels%B(2, ic) = B12(l) * factor

                    call apply_periodic_bc(parcels%position(:, ic))
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
                call parcel_replace(isma(m), l)

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
