module parcel_point
    use parcel_container, only : parcels, n_parcels
    use parcel_merge, only : pack_parcels
    use constants, only : zero, one, two , four, fpi, f12
    use parcel_interpl, only : trilinear, ngp
    use fields, only : velgradg, volg
    use parcel_bc
    use options, only : verbose
    use parameters, only : vcell, nx, nz, dx, lower, upper

    contains

        function get_eigenvalue(S) result(eval)
            double precision, intent(in) :: S(4)
            double precision             :: eval
            eval = dsqrt(dabs((S(1) - S(4)) ** 2 + four * S(2) * S(3)))
            eval = f12 * (S(1) + S(4) + eval)
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
            max_stretch = prefactor * dlog(threshold)

            do n = 1, last_index
                if (.not. parcels%stretch(n) > max_stretch) then
                    cycle
                endif

                ! evaluate velocity gradient at parcel locations
                S = get_velocity_gradient(parcels%position(n, :))

                ! separate along direction of eigenvector
                evec = get_eigenvector(S)

                parcels%volume(n) = f12 * parcels%volume(n)

                ! distance of child parcels to parent parcel
                h = parcels%stretch(n) * dsqrt(parcels%volume(n) * fpi)

                ! we only need to add one new parcel
                n_parcels = n_parcels + 1

                parcels%vorticity(n_parcels, :) = parcels%vorticity(n, :)
                parcels%volume(n_parcels) = parcels%volume(n)
                parcels%buoyancy(n_parcels) = parcels%buoyancy(n)
                parcels%humidity(n_parcels) = parcels%humidity(n)

                parcels%position(n_parcels, :) = parcels%position(n, :) - h * evec
                parcels%position(n, :) = parcels%position(n, :) + h * evec

                parcels%stretch(n) = zero
                parcels%stretch(n_parcels) = zero

                call apply_periodic_bc(parcels%position(n, :))
                call apply_periodic_bc(parcels%position(n_parcels, :))

                if (parcels%position(n, 2) > upper(2)) then
                    parcels%position(n, 2) = two * upper(2) &
                                           - parcels%position(n, 2)
                endif

                if (parcels%position(n_parcels, 2) < lower(2)) then
                    parcels%position(n_parcels, 2) = two * lower(2) &
                                                   - parcels%position(n_parcels, 2)
                endif
            enddo

            if (verbose) then
                print "(a36, i0, a3, i0)", &
                      "no. parcels before and after split: ", last_index, "...", n_parcels
            endif

        end subroutine point_split

        subroutine point_merge(vfraction)
            double precision, intent(in) :: vfraction
            double precision             :: vmin
            integer                      :: isma(0:n_parcels), n_remove, l
            double precision             :: vres(n_parcels), bres(n_parcels), zres(n_parcels)
            double precision             :: res_volg(-1:nz+1, 0:nx-1)
            double precision             :: ori_volg(-1:nz+1, 0:nx-1)
            double precision             :: res_vortg(-1:nz+1, 0:nx-1)
            double precision             :: res_tbuoyg(-1:nz+1, 0:nx-1)
            double precision             :: pos(2), weights(ngp), ww, volfi
            integer :: is(ngp), js(ngp)

            vmin = vcell / dble(vfraction)

            isma = zero

            res_volg = zero
            ori_volg = zero
            res_vortg = zero
            res_tbuoyg = zero

            n_remove = 0
            do n = 1, n_parcels
                pos = parcels%position(n, :)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call trilinear(pos, is, js, weights)

                if (parcels%volume(n) < vmin) then
                    n_remove = n_remove + 1
                    isma(n_remove) = n

                    do l = 1, ngp
                        ww = weights(l) * parcels%volume(n)
                        res_volg(js(l), is(l)) = res_volg(js(l), is(l)) + ww

                        res_vortg(js(l), is(l)) = res_vortg(js(l), is(l)) &
                                                + ww * parcels%vorticity(n, 1)

                        res_tbuoyg(js(l), is(l)) = res_tbuoyg(js(l), is(l)) &
                                                 + ww * parcels%buoyancy(n)
                    enddo
                else
                    do l = 1, ngp
                        ww = weights(l) * parcels%volume(n)
                        ori_volg(js(l), is(l)) = ori_volg(js(l), is(l)) + ww
                    enddo
                endif
            enddo


            ! edge-value doubling
            ori_volg(0,  :) = two * ori_volg(0,  :)
            ori_volg(nz, :) = two * ori_volg(nz, :)

            res_volg(0,  :) = two * res_volg(0,  :)
            res_volg(nz, :) = two * res_volg(nz, :)

            res_vortg(0,  :) = two * res_vortg(0,  :)
            res_vortg(nz, :) = two * res_vortg(nz, :)

            res_tbuoyg(0,  :) = two * res_tbuoyg(0,  :)
            res_tbuoyg(nz, :) = two * res_tbuoyg(nz, :)

            if (n_remove > 0) then
                ! remove invalid parcels
                call pack_parcels(isma, n_remove)

                do n = 1, n_parcels
                    vres(n) = zero
                    zres(n) = zero
                    bres(n) = zero

                    pos = parcels%position(n, :)
                    ! ensure parcel is within the domain
                    call apply_periodic_bc(pos)

                    ! get interpolation weights and mesh indices
                    call trilinear(pos, is, js, weights)

                    do l = 1, ngp
                        ww = weights(l) * parcels%volume(n) / ori_volg(js(l), is(l))
                        vres(n) = vres(n) &
                                + ww * res_volg(js(l), is(l))

                        bres(n) = bres(n) &
                                + ww * res_tbuoyg(js(l), is(l))

                        zres(n) = zres(n) &
                                + ww * res_vortg(js(l), is(l))
                    enddo
                    volfi = one / (parcels%volume(n) + vres(n))
                    parcels%buoyancy(n) = (parcels%buoyancy(n) * parcels%volume(n) + bres(n)) * volfi
                    parcels%vorticity(n, 1) = (parcels%vorticity(n, 1) * parcels%volume(n) + zres(n)) * volfi
                    parcels%volume(n) = parcels%volume(n) + vres(n)
                enddo
            endif
        end subroutine point_merge

end module parcel_point
