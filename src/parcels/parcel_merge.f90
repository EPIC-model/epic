module parcel_merge
    use nearest
    use constants, only : pi, max_num_parcels
    use parcel_container, only : parcel_container_type, n_parcels, parcel_replace
    use ellipse, only : get_B22
    use parameters, only : parcel_info
    implicit none

    private :: geometric_merge, &
               replace_parcels

    contains
        subroutine merge_ellipses(parcels)
            type(parcel_container_type), intent(inout) :: parcels
            integer                                    :: isma(max_num_parcels)
            integer                                    :: ibig(max_num_parcels)
            integer                                    :: n_merge ! number of merges

            ! find parcels to merge
            call find_nearest(isma, ibig, n_merge)

            if (n_merge > 0) then
                ! merge small parcels into large parcels
                if (parcel_info%merge_type == 'geometric') then
                    call geometric_merge(parcels, isma(1:n_merge), ibig(1:n_merge))
                else if (parcel_info%merge_type == 'optimal') then
                    call optimal_merge(parcels, isma(1:n_merge), ibig(1:n_merge))
                else
                    print *, "Unknown merge type '", trim(parcel_info%merge_type), "'."
                endif

                ! overwrite invalid parcels
                call replace_parcels(isma)
            endif

        end subroutine merge_ellipses

        ! merge ith parcel into jth parcel
        subroutine geometric_merge(parcels, isma, ibig)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(:)
            integer,                     intent(in)    :: ibig(:)
            integer                                    :: n, i, j
            double precision                           :: B11_1, B11_2
            double precision                           :: B12_1, B12_2
            double precision                           :: B22_1, B22_2, B22
            double precision                           :: a1b1, a2b2, ab, isqrab
            double precision                           :: mu1, mu2, zet, eta
            double precision                           :: mu11, mu22, mu12, detB

            do n = 1, size(isma)
                i = isma(n)
                j = ibig(n)

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

                parcels%B(j, 1) = mu12 * zet ** 2  + mu11 * B11_1 + mu22 * B11_2
                parcels%B(j, 2) = mu12 * zet * eta + mu11 * B12_1 + mu22 * B12_2
                B22             = mu12 * eta ** 2  + mu11 * B22_2 + mu22 * B22_2

                ! normalize such that determinant of the merger is 1
                detB = parcels%B(j, 1) * B22 - parcels%B(j, 2) ** 2

                parcels%B(j, 1) = parcels%B(j, 1) / detB
                parcels%B(j, 2) = parcels%B(j, 2) / detB
            enddo
        end subroutine geometric_merge


        subroutine optimal_merge(parcels, isma, ibig)
            type(parcel_container_type), intent(inout) :: parcels
            integer,                     intent(in)    :: isma(:)
            integer,                     intent(in)    :: ibig(:)


        end subroutine optimal_merge


        ! move invalid parcels to the end and update n_parcels
        subroutine replace_parcels(isma)
            integer, intent(in) :: isma(:)
            integer             :: l, n, m, k

            ! last entry we have to check
            k = size(isma)

            ! l points always to the last valid parcel
            l = n_parcels

            ! already update the number of parcels
            n_parcels = n_parcels - k

            ! find last parcel which is not invalid
            do while (l == isma(k) .and. k > 1)
                l = l - 1
                k = k - 1
            enddo

            if (l == 0) then
                print *, "Error: All parcels are invalid."
                stop
            endif

            ! move invalid parcels to the end of the container
            ! and update the n_parcels count
            m = 1
            n = 1

            do while (m <= k)
                if (n .ne. isma(m)) then
                    ! continue to next iteration if
                    ! this parcel is valid
                    n = n + 1
                    cycle
                endif

                ! invalid parcel; swap *this* with last valid parcel
                call parcel_replace(n, l)

                l = l - 1

                ! find next valid last parcel
                do while (l == isma(k) .and. k > 1)
                    l = l - 1
                    k = k - 1
                enddo

                ! next invalid
                m = m + 1

                n = n + 1
            enddo
        end subroutine replace_parcels

end module parcel_merge
