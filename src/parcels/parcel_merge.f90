module parcel_merge
    use nearest
    use constants, only : pi, max_num_parcels
    use parcel_container, only : parcel_container_type, n_parcels, parcel_replace
    use ellipse, only : get_B22
    use parameters, only : parcel_info, verbose
    implicit none

    private :: geometric_merge, &
               do_merge,        &
               optimal_merge,   &
               pack_parcels,    &
               remove_isolated

    contains
        subroutine merge_ellipses(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            integer                                    :: isma(max_num_parcels)
            integer                                    :: ibig(max_num_parcels)
            integer                                    :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(isma, ibig, n_merge)

            ! remove isolated parcels and update n_merge accordingly
            call remove_isolated(isma, ibig, n_merge)

            if (verbose) then
                print "(a36, i0, a3, i0)",                               &
                      "no. parcels before and after merge: ", n_parcels, &
                      "...", n_parcels - n_merge
            endif

            if (n_merge > 0) then
                ! merge small parcels into large parcels
                if (parcel_info%merge_type == 'geometric') then
                    call geometric_merge(parcels, isma, ibig, n_merge)
                else if (parcel_info%merge_type == 'optimal') then
                    call optimal_merge(parcels, isma, ibig, n_merge)
                else
                    print *, "Unknown merge type '", trim(parcel_info%merge_type), "'."
                endif

                ! overwrite invalid parcels
                call pack_parcels(isma, n_merge)
            endif

        end subroutine merge_ellipses

        ! remove isolated parcels;
        ! an isolated parcel is a parcel for which no big parcel is found (ibig = 0);
        subroutine remove_isolated(isma, ibig, n_merge)
            integer, intent(inout) :: isma(:), ibig(:)
            logical, allocatable   :: mask(:)
            integer, intent(out)   :: n_merge

            mask = ibig /= 0

            ! count number of .true. values
            n_merge = count(mask)

            ibig = pack(ibig, mask)
            isma = pack(isma, mask)

        end subroutine remove_isolated

        ! merge ith parcel into jth parcel (without B matrix scaling)
        subroutine do_merge(parcels, i, j, B11, B12, B22)
            type(parcel_container_type), intent(inout)  :: parcels
            integer,                     intent(in)     :: i, j
            double precision,            intent(out)    :: B11, B12, B22
            double precision                            :: B11_1, B11_2
            double precision                            :: B12_1, B12_2
            double precision                            :: B22_1, B22_2
            double precision                            :: a1b1, a2b2, isqrab, ab
            double precision                            :: mu1, mu2, zet, eta
            double precision                            :: mu11, mu22, mu12

            B11_1 = parcels%B(i, 1)
            B11_2 = parcels%B(j, 1)

            B12_1 = parcels%B(i, 2)
            B12_2 = parcels%B(j, 2)

            B22_1 = get_B22(B11_1, B12_1)
            B22_2 = get_B22(B11_2, B12_2)

            a1b1 = parcels%volume(i, 1) / pi
            a2b2 = parcels%volume(j, 1) / pi

            ab = a1b1 + a2b2
            isqrab = 1.0 / sqrt(ab)

            mu1 = a1b1 / ab
            mu2 = a2b2 / ab

            zet = 2.0 * isqrab * (parcels%position(j, 1) - parcels%position(i, 1))
            eta = 2.0 * isqrab * (parcels%position(j, 2) - parcels%position(i, 2))

            mu12 = mu1 * mu2
            mu11 = mu1 * mu1
            mu22 = mu2 * mu2

            B11 = mu12 * zet ** 2  + mu11 * B11_1 + mu22 * B11_2
            B12 = mu12 * zet * eta + mu11 * B12_1 + mu22 * B12_2
            B22 = mu12 * eta ** 2  + mu11 * B22_2 + mu22 * B22_2


            ! update center of mass
            parcels%position(j, 1) = mu1 * parcels%position(i, 1) &
                                   + mu2 * parcels%position(j, 1)

            parcels%position(j, 2) = mu1 * parcels%position(i, 2) &
                                   + mu2 * parcels%position(j, 2)

            ! update volume
            parcels%volume(j, 1) = ab * pi
        end subroutine do_merge


        subroutine geometric_merge(parcels, isma, ibig, n_merge)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(:)
            integer,                     intent(in)    :: ibig(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: n, i, j
            double precision                           :: B11, B12, B22, detB

            do n = 1, n_merge
                i = isma(n)
                j = ibig(n)

                ! merge small into big parcel --> return B11, B12, B22
                call do_merge(parcels, i, j, B11, B12, B22)

                ! normalize such that determinant of the merger is 1
                detB = B11 * B22 - B12 ** 2

                parcels%B(j, 1) = B11 / detB
                parcels%B(j, 2) = B12 / detB
            enddo
        end subroutine geometric_merge


        subroutine optimal_merge(parcels, isma, ibig, n_merge)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(:)
            integer,                     intent(in)    :: ibig(:)
            integer,                     intent(in)    :: n_merge
            integer                                    :: n, i, j
            double precision                           :: B11, B12, B22
            double precision                           :: mu, detB, merr, mup
            double precision                           :: a, b ,c

            do n = 1, n_merge
                i = isma(n)
                j = ibig(n)

                ! merge small into big parcel --> return B11, B12, B22
                call do_merge(parcels, i, j, B11, B12, B22)

                ! Solve the quartic to find best fit ellipse:
                !
                !
                !   Newton-Raphson to get smallest root:
                !      mu_{n+1} = mu_{n} - f(mu_{n}) / f'(mu_{n})
                !
                !      where
                !          f(mu_{n})  = mu_{n} ** 4 + b * mu_{n} ** 2 + a * mu_{n} + c
                !          f'(mu_{n}) = 4 * mu_{n} ** 3 + b * 2 * mu_{n} + a
                detB = B11 * B22 - B12 * B12
                a = B11 ** 2 + B22 ** 2 + 2.0 * B12 ** 2
                b = -2.0 - detB
                c = 1.0 - detB

                ! initial guess
                mu = - c / a

                merr = 1.0
                do while (merr > 1.e-12)
                    mup = (c + mu * (a + mu * (b + mu * mu))) / (a + mu * (2.0 * b + 4.0 * mu * mu))
                    mu = mu - mup
                    merr = abs(mup)
                enddo

                ! optimal B
                parcels%B(j, 1) = (B11 - mu * B22) / (1.0 - mu ** 2)
                parcels%B(j, 2) = B12 / (1.0 - mu)
            enddo
        end subroutine optimal_merge


        ! this algorithm replaces invalid parcels with valid parcels
        ! from the end of the container
        subroutine pack_parcels(isma, n_merge)
            integer, intent(in) :: isma(:)
            integer, intent(in) :: n_merge
            integer             :: k, l, n, m

            ! l points always to the last valid parcel in x
            l = n_parcels

            ! k points always to last invalid parcel in isma
            k = n_merge

            ! already update number of valid parcels
            n_parcels = n_parcels - k

            ! find last parcel which is not invalid
            do while (k > 0 .and. l == isma(k))
                l = l - 1
                k = k - 1
            enddo

            if (l == 0) then
                print *, "Error: All parcels are invalid."
                stop
            endif

            ! move invalid parcels to the end of the container
            m = 1

            ! n iterates over x array
            n = 1

            do while (m <= k)
                if (n .ne. isma(m)) then
                    ! continue to next iteration if
                    ! this parcel is valid
                    n = n + 1
                    cycle
                endif

                ! invalid parcel; overwrite *n* with last valid parcel *l*
                call parcel_replace(n, l)

                l = l - 1

                ! find next valid last parcel
                do while (k > 0 .and. l == isma(k))
                    l = l - 1
                    k = k - 1
                enddo


                ! next invalid
                m = m + 1

                n = n + 1
            enddo
        end subroutine pack_parcels


end module parcel_merge
