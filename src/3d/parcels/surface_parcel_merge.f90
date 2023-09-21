! =============================================================================
!                       Module to merge ellipses
!           The module implements the geometric merge procedure.
! =============================================================================
module surface_parcel_merging
    use surface_parcel_nearest
    use constants, only : pi, zero, one, two, four
    use surface_parcel_container, only : surface_parcel_container_type  &
                                       , n_lo_surf_parcels              &
                                       , n_up_surf_parcels              &
                                       , lo_surf_parcels                &
                                       , up_surf_parcels                &
                                       , surface_parcel_replace         &
                                       , get_delx                       &
                                       , get_dely
    use parcel_ellipse, only : get_ab
    use options, only : parcel, verbose
    use surface_parcel_bc
    use timer, only : start_timer, stop_timer
    implicit none

    integer:: surf_merge_timer

    private :: geometric_merge,     &
               do_group_merge,      &
               do_merge_ellipses,   &
               pack_parcels

    contains

        subroutine surface_parcel_merge
            call do_merge_ellipses(lo_surf_parcels, n_lo_surf_parcels)
            call do_merge_ellipses(up_surf_parcels, n_up_surf_parcels)
        end subroutine surface_parcel_merge

        ! Merge small parcels into neighbouring equal-sized parcels or bigger
        ! parcels which are close by.
        ! @param[inout] parcels is the parcel container
        subroutine do_merge_ellipses(s_parcels, n_par)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
            integer, allocatable, dimension(:)                 :: isma
            integer, allocatable, dimension(:)                 :: iclo
            integer                                            :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(s_parcels, n_par, isma, iclo, n_merge)

            call start_timer(surf_merge_timer)

            if (n_merge > 0) then
                ! merge small parcels into other parcels
                call geometric_merge(s_parcels, n_par, isma, iclo, n_merge)

                ! overwrite invalid parcels
                call pack_parcels(s_parcels, n_par, isma, n_merge)
            endif

            if (allocated(isma)) then
                deallocate(isma)
                deallocate(iclo)
            endif

            call stop_timer(surf_merge_timer)

        end subroutine do_merge_ellipses


        ! Actual merge.
        ! @param[inout] s_parcels is the surface parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        ! @param[out] B11m are the B11 matrix entries of the mergers
        ! @param[out] B12m are the B12 matrix entries of the mergers
        ! @param[out] B22m are the B22 matrix entries of the mergers
        ! @param[out] am are the areas of the mergers
        subroutine do_group_merge(s_parcels, n_par, isma, iclo, n_merge, B11m, B12m, B22m, am)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(in)    :: n_par
            integer,                             intent(in)    :: isma(0:)
            integer,                             intent(in)    :: iclo(:)
            integer,                             intent(in)    :: n_merge
            integer                                            :: m, ic, is, l, n
            integer                                            :: loca(n_par)
            double precision                                   :: pos0(2, n_merge)
            double precision                                   :: posm(2, n_merge), delx
            double precision                                   :: amerge, dely, mu
            double precision,                    intent(out)   :: B11m(n_merge), B12m(n_merge), &
                                                                  B22m(n_merge), am(n_merge)

            loca = zero

            l = 0
            do m = 1, n_merge
                ic = iclo(m) ! Index of closest other parcel

                if (loca(ic) == 0) then
                    ! Start a new merged parcel, indexed l:
                    l = l + 1
                    loca(ic) = l

                    ! am will contain the total area of the merged parcel
                    am(l) = s_parcels%area(ic)

                    !pos0 stores the x and y centre of the other parcel
                    pos0(:, l) = s_parcels%position(:, ic)

                    ! posm(:, l) will sum v(is)*(x(is)-x(ic)) and v(is)*(y(is)-y(ic))
                    ! modulo periodicity
                    posm(:, l) = zero

                    B11m(l) = zero
                    B12m(l) = zero
                    B22m(l) = zero
                endif

                ! Sum up all the small parcels merging with a common other one:
                ! "is" refers to the small parcel index
                is = isma(m) !Small parcel
                n = loca(ic)  !Index of merged parcel
                am(n) = am(n) + s_parcels%area(is) !Accumulate area of merged parcel

                ! works across periodic edge
                delx = get_delx(s_parcels%position(1, is), pos0(1, n))
                dely = get_dely(s_parcels%position(2, is), pos0(2, n))

                ! Accumulate sum of v(is)*(x(is)-x(ic)) and v(is)*(y(is)-y(ic))
                posm(1, n) = posm(1, n) + s_parcels%area(is) * delx
                posm(2, n) = posm(2, n) + s_parcels%area(is) * dely

            enddo

            ! Obtain the merged parcel centres
            ! (l = total number of merged parcels)
            do m = 1, l
                ! temporary scalar containing 1 / am(m)
                amerge = one / am(m)

                ! need to sanitise input and output, but first to determine input
                posm(:, m) = - amerge * posm(:, m)

                call apply_surface_periodic_bc(posm(:, m))

                ! x and y centre of merged parcel, modulo periodicity
                posm(1, m) = get_delx(pos0(1, m), posm(1, m))
                posm(2, m) = get_dely(pos0(2, m), posm(2, m))

                ! need to correct position
                call apply_surface_periodic_bc(posm(:, m))

            enddo

            loca = zero
            l = 0

            do m = 1, n_merge
                ic = iclo(m)

                if (loca(ic) == 0) then
                    l = l + 1
                    loca(ic) = l

                    amerge = one / am(l)

                    delx = get_delx(s_parcels%position(1, ic), posm(1, l))
                    dely = get_dely(s_parcels%position(2, ic), posm(2, l))

                    mu = s_parcels%area(ic) * amerge
                    B11m(l) = mu * (four * delx ** 2 + s_parcels%B(1, ic))
                    B12m(l) = mu * (four * delx * dely + s_parcels%B(2, ic))
                    B22m(l) = mu * (four * dely ** 2 + s_parcels%B(3, ic))

                    s_parcels%area(ic)  = am(l)
                    s_parcels%position(:, ic) = posm(:, l)

                endif

                is = isma(m)
                n = loca(ic)

                amerge = one / am(n)

                delx = get_delx(s_parcels%position(1, is), posm(1, n))
                dely = get_dely(s_parcels%position(2, is), posm(2, n))

                ! area fraction A_{is} / A
                mu = amerge * s_parcels%area(is)

                B11m(n) = B11m(n) + mu * (four * delx ** 2   + s_parcels%B(1, is))
                B12m(n) = B12m(n) + mu * (four * delx * dely + s_parcels%B(2, is))
                B22m(n) = B22m(n) + mu * (four * dely ** 2   + s_parcels%B(3, is))
            enddo

        end subroutine do_group_merge


        ! Geometric merging -- called by subroutine do_merge_ellipses.
        ! @param[inout] s_parcels is the surface parcel container
        ! @param[in] isma are the indices of the small parcels
        ! @param[in] iclo are the indices of the close parcels
        ! @param[in] n_merge is the array size of isma and iclo
        subroutine geometric_merge(s_parcels, n_par, isma, iclo, n_merge)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                            intent(in)     :: n_par
            integer,                            intent(in)     :: isma(0:)
            integer,                            intent(in)     :: iclo(:)
            integer,                            intent(in)     :: n_merge
            integer                                            :: m, ic, l
            integer                                            :: loca(n_par)
            double precision                                   :: factor
            double precision                                   :: B11(n_merge), &
                                                                  B12(n_merge), &
                                                                  B22(n_merge), &
                                                                  A(n_merge)

            call do_group_merge(s_parcels, n_par, isma, iclo, n_merge, B11, B12, B22, A)

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
                    factor = get_ab(A(l)) / dsqrt(B11(l) * B22(l) - B12(l) ** 2)

                    s_parcels%B(1, ic) = B11(l) * factor
                    s_parcels%B(2, ic) = B12(l) * factor
                    s_parcels%B(3, ic) = B22(l) * factor

                    call apply_surface_periodic_bc(s_parcels%position(:, ic))
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
        subroutine pack_parcels(s_parcels, n_par, isma, n_merge)
            type(surface_parcel_container_type), intent(inout) :: s_parcels
            integer,                             intent(inout) :: n_par
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
                call surface_parcel_replace(s_parcels, isma(m), l)

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


end module surface_parcel_merging
