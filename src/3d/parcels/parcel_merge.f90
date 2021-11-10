! =============================================================================
!                       Module to merge ellipsoids
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_merge
    use parcel_nearest
    use constants, only : pi, zero, one, two, five, f13 &
                        , max_num_parcels
    use parcel_container, only : parcel_container_type  &
                               , n_parcels              &
                               , parcel_replace         &
                               , get_delx               &
                               , get_dely
    use parcel_ellipsoid, only : get_B33, get_abc
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
        ! @param[out] Bm are the B matrix entries of the mergers
        ! @param[out] vm are the volumes of the mergers
        subroutine do_group_merge(parcels, isma, iclo, n_merge, Bm, vm)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, is, l, n
            integer                                    :: loca(n_parcels)
            double precision                           :: x0(n_merge), y0(n_merge)
            double precision                           :: xm(n_merge), ym(n_merge), zm(n_merge)
            double precision                           :: delx, vmerge, dely, delz, B33, mu
            double precision                           :: buoym(n_merge), vortm(n_merge, 3)
#ifndef ENABLE_DRY_MODE
            double precision                           :: hum(n_merge)
#endif
            double precision,            intent(out)   :: Bm(n_merge, 6) ! B11, B12, B13, B22, B23, B33
            double precision,            intent(out)   :: vm(n_merge)

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
                    x0(l) = parcels%position(ic, 1)

                    !y0 stores the y centre of the other parcel
                    y0(l) = parcels%position(ic, 1)

                    ! xm will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    xm(l) = zero

                    ! ym will sum v(is)*(y(is)-y(ic)) modulo periodicity
                    ym(l) = zero

                    ! zm will contain v(ic)*z(ic)+sum{v(is)*z(is)}
                    zm(l) = parcels%volume(ic) * parcels%position(ic, 3)

                    ! buoyancy and humidity
                    buoym(l) = parcels%volume(ic) * parcels%buoyancy(ic)
#ifndef ENABLE_DRY_MODE
                    hum(l) = parcels%volume(ic) * parcels%humidity(ic)
#endif
                    vortm(l, :) = parcels%volume(ic) * parcels%vorticity(ic, :)

                    Bm(l, :) = zero
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                vm(n) = vm(n) + parcels%volume(is) !Accumulate volume of merged parcel

                ! works across periodic edge
                delx = get_delx(parcels%position(is, 1), x0(n))
                dely = get_dely(parcels%position(is, 2), y0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic)) and v(is)*(y(is)-y(ic))
                xm(n) = xm(n) + parcels%volume(is) * delx
                ym(n) = ym(n) + parcels%volume(is) * dely

                ! Accumulate v(ic)*z(ic)+sum{v(is)*z(is)}
                zm(n) = zm(n) + parcels%volume(is) * parcels%position(is, 3)

                ! Accumulate buoyancy and humidity
                buoym(n) = buoym(n) + parcels%volume(is) * parcels%buoyancy(is)
#ifndef ENABLE_DRY_MODE
                hum(n) = hum(n) + parcels%volume(is) * parcels%humidity(is)
#endif
                vortm(n, :) = vortm(n, :) + parcels%volume(is) * parcels%vorticity(is, :)
            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / vm(m)
                vmerge = one / vm(m)

                ! x and y centre of merged parcel, modulo periodicity
                xm(m) = get_delx(x0(m), - vmerge * xm(m))
                ym(m) = get_dely(y0(m), - vmerge * ym(m))

                ! z centre of merged parcel
                zm(m) = vmerge * zm(m)

                ! buoyancy and humidity
                buoym(m) = vmerge * buoym(m)
#ifndef ENABLE_DRY_MODE
                hum(m) = vmerge * hum(m)
#endif
                vortm(m, :) = vmerge * vortm(m, :)
            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    vmerge = one / vm(l)

                    B33 = get_B33(parcels%B(ic, 1), parcels%volume(ic))

                    delx = get_delx(parcels%position(ic, 1), xm(l))
                    dely = get_delx(parcels%position(ic, 2), ym(l))
                    delz = parcels%position(ic, 3) - zm(l)

                    mu = parcels%volume(ic) * vmerge

                    Bm(l, 1) = mu * (five * delx ** 2   + parcels%B(ic, 1))
                    Bm(l, 2) = mu * (five * delx * dely + parcels%B(ic, 2))
                    Bm(l, 3) = mu * (five * delx * delz + parcels%B(ic, 3))
                    Bm(l, 4) = mu * (five * dely ** 2   + parcels%B(ic, 4))
                    Bm(l, 5) = mu * (five * dely * delz + parcels%B(ic, 5))
                    Bm(l, 6) = mu * (five * delz ** 2   + B33)

                    parcels%volume(ic)  = vm(l)
                    parcels%position(ic, 1) = xm(l)
                    parcels%position(ic, 2) = ym(l)
                    parcels%position(ic, 3) = zm(l)

                    parcels%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ic) = hum(l)
#endif
                    parcels%vorticity(ic, :) = vortm(l, :)

                endif

                is = isma(m)
                n = loca(ic)

                vmerge = one / vm(n)

                delx = get_delx(parcels%position(is, 1), xm(n))
                dely = get_dely(parcels%position(is, 2), ym(n))
                delz = parcels%position(is, 3) - zm(n)

                B33 = get_B33(parcels%B(is, :), parcels%volume(is))

                ! volume fraction A_{is} / A
                mu = vmerge * parcels%volume(is)

                Bm(n, 1) = Bm(n, 1) + mu * (five * delx ** 2   + parcels%B(is, 1))
                Bm(n, 2) = Bm(n, 2) + mu * (five * delx * dely + parcels%B(is, 2))
                Bm(n, 3) = Bm(n, 3) + mu * (five * delx * delz + parcels%B(is, 3))
                Bm(n, 4) = Bm(n, 4) + mu * (five * dely ** 2   + parcels%B(is, 4))
                Bm(n, 5) = Bm(n, 5) + mu * (five * dely * delz + parcels%B(is, 5))
                Bm(n, 6) = Bm(n, 6) + mu * (five * delz ** 2   + B33)
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
            double precision                           :: factor, detB
            double precision                           :: B(n_merge, 6), &
                                                          V(n_merge)

            call do_group_merge(parcels, isma, iclo, n_merge, B, V)

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
                    detB = B(l, 1) * (B(l, 2) * B(l, 6) - B(l, 5) ** 2) &
                         - B(l, 2) * (B(l, 2) * B(l, 6) - B(l, 3) * B(l, 5)) &
                         + B(l, 3) * (B(l, 2) * B(l, 5) - B(l, 3) * B(l, 4))

                    factor = (get_abc(V(l)) ** 2 / detB) ** f13

                    parcels%B(ic, :) = B(l, 1:5) * factor

                    call apply_periodic_bc(parcels%position(ic, :))
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
