module parcel_point
    use parcel_container, only : parcels, n_parcels
    use constants, only : zero, one, two , four, fpi
    use parcel_interpl, only : trilinear, ngp
    use fields, only : velgradg
    use parcel_bc
    use options, only : verbose

    contains

        function get_eigenvalue(S) result(eval)
            double precision, intent(in) :: S(4)
            double precision             :: eval
            eval = dsqrt(dabs(S(1) ** 2 - two * S(1) * S(4) + four * S(2) * S(3) + S(4) ** 2))
            eval = 0.5d0 * (S(1) + S(4) + eval)
        end function get_eigenvalue

        function get_eigenvector(S) result(evec)
            double precision, intent(in) :: S(4)
            double precision             :: evec(2)
            double precision             :: eval

            eval = get_eigenvalue(S)
            evec(1) = (eval - S(4)) / S(3)
            evec(2) = one

            if (evec(1) == zero) then
                evec(1) = epsilon(evec(1))
            endif

            evec = evec / norm2(evec)
        end function get_eigenvector

        function get_max_stretch(threshold, prefactor) result(max_stretch)
            double precision, intent(in) :: threshold, prefactor
            double precision             :: max_stretch
            max_stretch = prefactor * dlog(threshold)
        end function get_max_stretch

        function get_velocity_gradient(inpos) result(S)
            double precision, intent(in) :: inpos(2)
            double precision             :: pos(2)
            double precision             :: S(4)
            double precision             :: weights(ngp)
            integer                      :: c, l, is(ngp), js(ngp)


            pos = inpos

            ! ensure parcel is within the domain
            call apply_periodic_bc(pos)

            call trilinear(pos, is, js, weights)

            S = zero

            do l = 1, ngp
                do c = 1, 4
                    S(c) = S(c) + weights(l) * velgradg(js(l), is(l), c)
                enddo
            enddo
        end function get_velocity_gradient

        ! Splitting according to the MPIC paper
        ! https://doi.org/10.1002/qj.3319
        subroutine point_split(threshold, prefactor)
            double precision, intent(in) :: threshold, prefactor
            double precision             :: max_stretch, evec(2), S(4)
            double precision             :: h
            integer                      :: n, last_index

            last_index = n_parcels
            max_stretch = get_max_stretch(threshold, prefactor)



            do n = 1, last_index
                if (.not. parcels%stretch(n) > max_stretch) then
                    cycle
                endif

                ! evaluate velocity gradient at parcel locations
                S = get_velocity_gradient(parcels%position(n, :))

                ! separate along direction of eigenvector
                evec = get_eigenvector(S)

                h = max_stretch * dsqrt(0.5d0 * parcels%volume(n) * fpi)

                parcels%volume(n) = 0.5d0 * parcels%volume(n)

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1

                parcels%velocity(n_parcels, :) = parcels%velocity(n, :)
                parcels%vorticity(n_parcels, :) = parcels%vorticity(n, :)
                parcels%volume(n_parcels) = parcels%volume(n)
                parcels%buoyancy(n_parcels) = parcels%buoyancy(n)
                parcels%humidity(n_parcels) = parcels%humidity(n)

                parcels%position(n_parcels, :) = parcels%position(n, :) - h * evec
                parcels%position(n, :) = parcels%position(n, :)  + h * evec

                parcels%stretch(n) = zero
                parcels%stretch(n_parcels) = zero

            enddo

            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", last_index, "...", n_parcels
            endif

        end subroutine point_split

end module parcel_point
