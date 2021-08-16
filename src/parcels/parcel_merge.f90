! =============================================================================
!                       Module to merge ellipses
! The module provides a 'geometric' and 'optimal' merge procedure which can be
! selected by the user at runtime.
! =============================================================================
module parcel_merge
    use parcel_nearest
    use constants, only : pi, zero, one, two, four    &
                        , max_num_parcels
    use parcel_container, only : parcel_container_type  &
                               , n_parcels              &
                               , parcel_replace         &
                               , get_delx
    use parcel_ellipse, only : get_B22, get_ab
    use options, only : parcel, verbose
    use parcel_bc
    use timer, only : start_timer, stop_timer

    implicit none

    integer :: merge_timer

    private :: geometric_merge, &
               optimal_merge,   &
               do_group_merge,  &
               solve_quartic,   &
               pack_parcels

    contains
        subroutine merge_ellipses(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            integer, allocatable, dimension(:)         :: isma
            integer, allocatable, dimension(:)         :: ibig
            integer                                    :: n_merge ! number of merges

            call start_timer(merge_timer)

            ! find parcels to merge
            call find_nearest(isma, ibig, n_merge)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", n_parcels, &
                      "...", n_parcels - n_merge
            endif
#endif

            if (n_merge > 0) then
                ! merge small parcels into large parcels
                if (parcel%merge_type == 'geometric') then
                    call geometric_merge(parcels, isma, ibig, n_merge)
                else if (parcel%merge_type == 'optimal') then
                    call optimal_merge(parcels, isma, ibig, n_merge)
                else
                    print *, "Unknown merge type '", trim(parcel%merge_type), "'."
                    stop
                endif

                ! overwrite invalid parcels
                call pack_parcels(isma, n_merge)
            endif

            deallocate(isma)
            deallocate(ibig)

            call stop_timer(merge_timer)

        end subroutine merge_ellipses

        ! merge is-th parcel into ib-th parcel (without B matrix scaling)
        subroutine do_bimerge(parcels, is, ib, B11, B12, B22, ab)
            type(parcel_container_type), intent(inout)  :: parcels
            integer,                     intent(in)     :: is, ib
            double precision,            intent(out)    :: B11, B12, B22, ab
            double precision                            :: B11_1, B11_2
            double precision                            :: B12_1, B12_2
            double precision                            :: B22_1, B22_2
            double precision                            :: a1b1, a2b2, abi
            double precision                            :: mu1, mu2, zet, eta, mu12
            double precision                            :: delx

            B11_1 = parcels%B(is, 1)
            B11_2 = parcels%B(ib, 1)

            B12_1 = parcels%B(is, 2)
            B12_2 = parcels%B(ib, 2)

            B22_1 = get_B22(B11_1, B12_1, parcels%volume(is))
            B22_2 = get_B22(B11_2, B12_2, parcels%volume(ib))

            a1b1 = get_ab(parcels%volume(is))
            a2b2 = get_ab(parcels%volume(ib))

            ab = a1b1 + a2b2
            abi = one / ab
            mu1 = a1b1 * abi
            mu2 = a2b2 * abi

            ! works across periodic edge
            delx = get_delx(parcels%position(ib, 1), parcels%position(is, 1))

            zet = two * delx

            eta = two * (parcels%position(ib, 2) - parcels%position(is, 2))

            mu12 = mu1 * mu2

            B11 = mu12 * zet ** 2  + mu1 * B11_1 + mu2 * B11_2
            B12 = mu12 * zet * eta + mu1 * B12_1 + mu2 * B12_2
            B22 = mu12 * eta ** 2  + mu1 * B22_1 + mu2 * B22_2

            ! update center of mass
            parcels%position(ib, 1) = - mu1 * delx &
                                    + parcels%position(ib, 1)

            parcels%position(ib, 2) = mu1 * parcels%position(is, 2) &
                                    + mu2 * parcels%position(ib, 2)

            ! update buoyancy, humidity and vorticity
            parcels%buoyancy(ib) = mu1 * parcels%buoyancy(is) &
                                 + mu2 * parcels%buoyancy(ib)
#ifndef ENABLE_DRY_MODE
            parcels%humidity(ib) = mu1 * parcels%humidity(is) &
                                 + mu2 * parcels%humidity(ib)
#endif
            parcels%vorticity(ib) = mu1 * parcels%vorticity(is) &
                                  + mu2 * parcels%vorticity(ib)

            ! update volume
            parcels%volume(ib) = ab * pi
        end subroutine do_bimerge


        subroutine do_group_merge(parcels, isma, ibig, n_merge, B11m, B12m, B22m, vm)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: ibig(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ib, is, l, n
            integer                                    :: loca(n_parcels)
            double precision                           :: x0(n_merge), xm(n_merge)
            double precision                           :: zm(n_merge), delx, vmerge, dely, B22, mu
            double precision                           :: buoym(n_merge), vortm(n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                           :: hum(n_merge)
#endif
            double precision,            intent(out)   :: B11m(n_merge), B12m(n_merge), B22m(n_merge), &
                                                          vm(n_merge)

            loca = zero

            l = 0
            do m = 1, n_merge
                ib = ibig(m) ! Index of large parcel

                if (loca(ib) .eq. 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ib) = l

                    ! vm will contain the total volume of the merged parcel
                    vm(l) = parcels%volume(ib)

                    !x0 stores the x centre of the large parcel
                    x0(l) = parcels%position(ib, 1)

                    ! xm will sum v(is)*(x(is)-x(ib)) modulo periodicity
                    xm(l) = zero

                    ! zm will contain v(ib)*z(ib)+sum{v(is)*z(is)}
                    zm(l) = parcels%volume(ib) * parcels%position(ib, 2)

                    ! buoyancy and humidity
                    buoym(l) = parcels%volume(ib) * parcels%buoyancy(ib)
#ifndef ENABLE_DRY_MODE
                    hum(l) = parcels%volume(ib) * parcels%humidity(ib)
#endif
                    vortm(l) = parcels%volume(ib) * parcels%vorticity(ib)

                    B11m(l) = zero
                    B12m(l) = zero
                    B22m(l) = zero
                endif

                ! Sum up all the small parcels merging with a common larger one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ib)  !Index of merged parcel
                vm(n) = vm(n) + parcels%volume(is) !Accumulate volume of merged parcel

                ! works across periodic edge
                delx = get_delx(parcels%position(is, 1), x0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ib))
                xm(n) = xm(n) + parcels%volume(is) * delx

                ! Accumulate v(ib)*z(ib)+sum{v(is)*z(is)}
                zm(n) = zm(n) + parcels%volume(is) * parcels%position(is, 2)

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
                xm(m) = get_delx(x0(m), - vmerge * xm(m))

                ! z centre of merged parcel
                zm(m) = vmerge * zm(m)

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
                ib = ibig(m)

                if (loca(ib) .eq. 0) then
                    l = l + 1
                    loca(ib) = l

                    vmerge = one / vm(l)

                    B22 = get_B22(parcels%B(ib, 1), parcels%B(ib, 2), parcels%volume(ib))

                    delx = get_delx(parcels%position(ib, 1), xm(l))
                    dely = parcels%position(ib, 2) - zm(l)

                    mu = parcels%volume(ib) * vmerge
                    B11m(l) = mu * (four * delx ** 2 + parcels%B(ib, 1))
                    B12m(l) = mu * (four * delx * dely + parcels%B(ib, 2))
                    B22m(l) = mu * (four * dely ** 2 + B22)

                    parcels%volume(ib)  = vm(l)
                    parcels%position(ib, 1) = xm(l)
                    parcels%position(ib, 2) = zm(l)

                    parcels%buoyancy(ib) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ib) = hum(l)
#endif
                    parcels%vorticity(ib) = vortm(l)

                endif

                is = isma(m)
                n = loca(ib)

                vmerge = one / vm(n)

                delx = get_delx(parcels%position(is, 1), xm(n))
                dely = parcels%position(is, 2) - zm(n)

                B22 = get_B22(parcels%B(is, 1), parcels%B(is, 2), parcels%volume(is))

                ! volume fraction A_{is} / A
                mu = vmerge * parcels%volume(is)

                B11m(n) = B11m(n) + mu * (four * delx ** 2   + parcels%B(is, 1))
                B12m(n) = B12m(n) + mu * (four * delx * dely + parcels%B(is, 2))
                B22m(n) = B22m(n) + mu * (four * dely ** 2   + B22)
            enddo

        end subroutine do_group_merge


        subroutine geometric_merge(parcels, isma, ibig, n_merge)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: ibig(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ib, l
            integer                                    :: loca(n_parcels)
            double precision                           :: factor
            double precision                           :: B11(n_merge), &
                                                          B12(n_merge), &
                                                          B22(n_merge), &
                                                          V(n_merge)

            call do_group_merge(parcels, isma, ibig, n_merge, B11, B12, B22, V)


            loca = zero

            l = 0
            do m = 1, n_merge
                ib = ibig(m)

                if (loca(ib) .eq. 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ib) = l

                    ! normalize such that determinant of the merger is (ab)**2
                    ! ab / sqrt(det(B))
                    factor = get_ab(V(l)) / dsqrt(B11(l) * B22(l) - B12(l) ** 2)

                    parcels%B(ib, 1) = B11(l) * factor
                    parcels%B(ib, 2) = B12(l) * factor

                    call apply_periodic_bc(parcels%position(ib, :))
                endif
            enddo

        end subroutine geometric_merge


        subroutine optimal_merge(parcels, isma, ibig, n_merge)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: ibig(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ib, l
            integer                                    :: loca(n_parcels)
            double precision                           :: mu, ab
            double precision                           :: B11(n_merge), &
                                                          B12(n_merge), &
                                                          B22(n_merge), &
                                                          V(n_merge)

            call do_group_merge(parcels, isma, ibig, n_merge, B11, B12, B22, V)

            loca = zero

            l = 0
            do m = 1, n_merge
                ib = ibig(m)

                if (loca(ib) .eq. 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ib) = l

                    ab = get_ab(V(l))
                    mu = solve_quartic(B11(l), B12(l), B22(l), ab)

                    ! optimal B
                    parcels%B(ib, 1) = (B11(l) - mu * B22(l)) / (one - mu ** 2)
                    parcels%B(ib, 2) = B12(l) / (one - mu)

                    call apply_periodic_bc(parcels%position(ib, :))
                endif
            enddo
        end subroutine optimal_merge


        function solve_quartic(B11, B12, B22, ab) result(mu)
            double precision, intent(in) :: B11, B12, B22, ab
            double precision             :: mu, detB, merr, mup
            double precision             :: a, b ,c, a2b2i
            ! Solve the quartic to find best fit ellipse:
            !
            !
            !   Newton-Raphson to get smallest root:
            !      mu_{n+1} = mu_{n} - f(mu_{n}) / f'(mu_{n})
            !
            !      where
            !          f(mu_{n})  = mu_{n} ** 4 + b * mu_{n} ** 2 + a * mu_{n} + c
            !          f'(mu_{n}) = 4 * mu_{n} ** 3 + b * 2 * mu_{n} + a
            a2b2i = one / ab ** 2
            detB = (B11 * B22 - B12 * B12) * a2b2i
            a = (B11 ** 2 + B22 ** 2 + two * B12 ** 2) * a2b2i
            b = -two - detB
            c = one - detB

            ! initial guess
            mu = - c / a

            merr = 1.0
            do while (merr > 1.e-12)
                mup = (c + mu * (a + mu * (b + mu * mu))) / (a + mu * (two * b + four * mu * mu))
                mu = mu - mup
                merr = abs(mup)
            enddo

        end function solve_quartic


        ! this algorithm replaces invalid parcels with valid parcels
        ! from the end of the container
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
