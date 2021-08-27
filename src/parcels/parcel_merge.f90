! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
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
#ifdef ENABLE_VERBOSE
    use merge_hdf5, only : write_h5_mergees,            &
                           write_h5_mergers,            &
                           write_h5_parcels_in_cell
#endif
    implicit none

    integer :: merge_timer

    private :: geometric_merge, &
               do_group_merge,  &
               pack_parcels

    contains
        subroutine merge_ellipses(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            integer, allocatable, dimension(:)         :: isma
            integer, allocatable, dimension(:)         :: iclo
            integer                                    :: n_merge ! number of merges

            call start_timer(merge_timer)

            ! find parcels to merge
            call find_nearest(isma, iclo, n_merge)

#ifdef ENABLE_VERBOSE
            if (verbose) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", n_parcels, &
                      "...", n_parcels - n_merge
            endif
#endif

            if (n_merge > 0) then
#ifdef ENABLE_MERGER_DUMP
                call write_h5_mergees(isma, iclo, n_merge)
                call write_h5_parcels_in_cell(iclo, n_merge)
#endif
                ! merge small parcels into other parcels
                call geometric_merge(parcels, isma, iclo, n_merge)

#ifdef ENABLE_MERGER_DUMP
                call write_h5_mergers(iclo, n_merge)
#endif

                ! overwrite invalid parcels
                call pack_parcels(isma, n_merge)
            endif

            if (allocated(isma)) then
                deallocate(isma)
                deallocate(iclo)
            endif

            call stop_timer(merge_timer)

        end subroutine merge_ellipses


        subroutine do_group_merge(parcels, isma, iclo, n_merge, B11m, B12m, B22m, vm)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(0:)
            integer,                     intent(in)    :: iclo(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: m, ic, is, l, n
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
                ic = iclo(m) ! Index of closest other parcel

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! vm will contain the total volume of the merged parcel
                    vm(l) = parcels%volume(ic)

                    !x0 stores the x centre of the other parcel
                    x0(l) = parcels%position(ic, 1)

                    ! xm will sum v(is)*(x(is)-x(ic)) modulo periodicity
                    xm(l) = zero

                    ! zm will contain v(ic)*z(ic)+sum{v(is)*z(is)}
                    zm(l) = parcels%volume(ic) * parcels%position(ic, 2)

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
                delx = get_delx(parcels%position(is, 1), x0(n))

                ! Accumulate sum of v(is)*(x(is)-x(ic))
                xm(n) = xm(n) + parcels%volume(is) * delx

                ! Accumulate v(ic)*z(ic)+sum{v(is)*z(is)}
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
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    vmerge = one / vm(l)

                    B22 = get_B22(parcels%B(ic, 1), parcels%B(ic, 2), parcels%volume(ic))

                    delx = get_delx(parcels%position(ic, 1), xm(l))
                    dely = parcels%position(ic, 2) - zm(l)

                    mu = parcels%volume(ic) * vmerge
                    B11m(l) = mu * (four * delx ** 2 + parcels%B(ic, 1))
                    B12m(l) = mu * (four * delx * dely + parcels%B(ic, 2))
                    B22m(l) = mu * (four * dely ** 2 + B22)

                    parcels%volume(ic)  = vm(l)
                    parcels%position(ic, 1) = xm(l)
                    parcels%position(ic, 2) = zm(l)

                    parcels%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ic) = hum(l)
#endif
                    parcels%vorticity(ic) = vortm(l)

                endif

                is = isma(m)
                n = loca(ic)

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
                    factor = get_ab(V(l)) / dsqrt(B11(l) * B22(l) - B12(l) ** 2)

                    parcels%B(ic, 1) = B11(l) * factor
                    parcels%B(ic, 2) = B12(l) * factor

                    call apply_periodic_bc(parcels%position(ic, :))
                endif
            enddo

        end subroutine geometric_merge


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
