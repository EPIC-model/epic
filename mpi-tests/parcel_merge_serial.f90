! =============================================================================
!                       Module to merge ellipsoids
!           The module implements the geometric merge procedure.
! =============================================================================
module parcel_merge_serial
    use parcel_nearest_serial
    use constants, only : pi, zero, one, two, five, f13
    use parcel_container, only : parcel_container_type      &
                               , n_parcels                  &
                               , parcel_replace             &
                               , get_delx_across_periodic   &
                               , get_dely_across_periodic
    use parcel_ellipsoid, only : get_B33, get_abc
    use options, only : parcel, verbose
    use parcel_bc
    use mpi_timer, only : start_timer, stop_timer

    implicit none

    integer :: merge_timer

    ! number of parcel merges (is reset in every write step)
    integer :: n_parcel_merges = 0

    private :: geometric_merge, &
               do_group_merge,  &
               pack_parcels

    contains

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        ! @param[inout] parcels is the parcel container
        subroutine merge_parcels(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            integer, allocatable, dimension(:)         :: isma
            integer, allocatable, dimension(:)         :: iclo
            integer                                    :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(isma, iclo, n_merge)

            n_parcel_merges = n_parcel_merges + n_merge

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

        end subroutine merge_parcels


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
            double precision                           :: posm(3, n_merge)
            double precision                           :: delx, vmerge, dely, delz, B33, mu
            double precision                           :: buoym(n_merge), vortm(3, n_merge)
#ifndef ENABLE_DRY_MODE
            double precision                           :: hum(n_merge)
#endif
            double precision,            intent(out)   :: Bm(6, n_merge) ! B11, B12, B13, B22, B23, B33
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
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
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

                    B33 = get_B33(parcels%B(:, ic), parcels%volume(ic))

                    delx = get_delx_across_periodic(parcels%position(1, ic), posm(1, l))
                    dely = get_dely_across_periodic(parcels%position(2, ic), posm(2, l))
                    delz = parcels%position(3, ic) - posm(3, l)

                    mu = parcels%volume(ic) * vmerge


                    Bm(1, l) = mu * (five * delx ** 2   + parcels%B(1, ic))
                    Bm(2, l) = mu * (five * delx * dely + parcels%B(2, ic))
                    Bm(3, l) = mu * (five * delx * delz + parcels%B(3, ic))
                    Bm(4, l) = mu * (five * dely ** 2   + parcels%B(4, ic))
                    Bm(5, l) = mu * (five * dely * delz + parcels%B(5, ic))
                    Bm(6, l) = mu * (five * delz ** 2   + B33)

                    parcels%volume(ic)  = vm(l)
                    parcels%position(1, ic) = posm(1, l)
                    parcels%position(2, ic) = posm(2, l)
                    parcels%position(3, ic) = posm(3, l)

                    parcels%buoyancy(ic) = buoym(l)
#ifndef ENABLE_DRY_MODE
                    parcels%humidity(ic) = hum(l)
#endif
                    parcels%vorticity(:, ic) = vortm(:, l)

                endif

                is = isma(m)
                n = loca(ic)

                vmerge = one / vm(n)

                delx = get_delx_across_periodic(parcels%position(1, is), posm(1, n))
                dely = get_dely_across_periodic(parcels%position(2, is), posm(2, n))
                delz = parcels%position(3, is) - posm(3, n)

                B33 = get_B33(parcels%B(:, is), parcels%volume(is))

                ! volume fraction V_{is} / V
                mu = vmerge * parcels%volume(is)

                Bm(1, n) = Bm(1, n) + mu * (five * delx ** 2   + parcels%B(1, is))
                Bm(2, n) = Bm(2, n) + mu * (five * delx * dely + parcels%B(2, is))
                Bm(3, n) = Bm(3, n) + mu * (five * delx * delz + parcels%B(3, is))
                Bm(4, n) = Bm(4, n) + mu * (five * dely ** 2   + parcels%B(4, is))
                Bm(5, n) = Bm(5, n) + mu * (five * dely * delz + parcels%B(5, is))
                Bm(6, n) = Bm(6, n) + mu * (five * delz ** 2   + B33)
            enddo

        end subroutine do_group_merge


        ! Geometric merging -- called by subroutine merge_parcels.
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
            double precision                           :: B(6, n_merge), &
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
                    detB = B(1, l) * (B(4, l) * B(6, l) - B(5, l) ** 2) &
                         - B(2, l) * (B(2, l) * B(6, l) - B(3, l) * B(5, l)) &
                         + B(3, l) * (B(2, l) * B(5, l) - B(3, l) * B(4, l))

                    factor = (get_abc(V(l)) ** 2 / detB) ** f13

                    parcels%B(:, ic) = B(1:5, l) * factor

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


end module parcel_merge_serial
