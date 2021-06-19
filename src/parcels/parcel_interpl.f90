! =============================================================================
! This module contains the subroutines to do parcel-to-grid and grid-to-parcel
! interpolation.
! =============================================================================
module parcel_interpl
    use constants, only : max_num_parcels, zero, one, two
    use parameters, only : nx, nz
    use options, only : parcel
    use parcel_container, only : parcels, n_parcels
    use parcel_bc, only : apply_periodic_bc
    use parcel_ellipse
    use fields
    implicit none

    private :: par2grid_elliptic,       &
               par2grid_non_elliptic,   &
               grid2par_elliptic,       &
               grid2par_non_elliptic


    ! number of indices and weights
    integer, parameter :: ngp = 4

    ! interpolation indices
    ! (first dimension x, y; second dimension k-th index)
    integer :: is(ngp), js(ngp)

    ! interpolation weights
    double precision :: weights(ngp)

    private :: is, js, weights

    contains

        subroutine vol2grid
            volg = zero

            if (parcel%is_elliptic) then
                call vol2grid_elliptic
            else
                call vol2grid_non_elliptic
            endif

            ! apply free slip boundary condition
            volg(0,  :) = two * volg(0,  :)
            volg(nz, :) = two * volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :) = volg(1,    :) + volg(-1,   :)
            volg(nz-1, :) = volg(nz-1, :) + volg(nz+1, :)

        end subroutine


        subroutine vol2grid_elliptic
            double precision  :: points(2, 2)
            integer           :: n, p, l
            double precision  :: pvol, pvor

            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(n, :), &
                                            pvol, parcels%B(n, :))


                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(p, :), is, js, weights)

                    do l = 1, ngp
                        volg(js(l), is(l)) = volg(js(l), is(l)) &
                                           + f12 * weights(l) * pvol
                    enddo
                enddo
            enddo
        end subroutine vol2grid_elliptic

        subroutine vol2grid_elliptic_symmetry_check
            double precision :: points(2, 2), V, B(2), pos(2)
            integer          :: n, p, l, m
            double precision :: pvol

            volg = zero

            do m = -1, 1, 2
                do n = 1, n_parcels

                    pos = parcels%position(n, :)
                    pvol = parcels%volume(n)
                    pos(1) = dble(m) * pos(1)
                    V = dble(m) * pvol
                    B = parcels%B(n, :)

                    B(2) = dble(m) * B(2)

                    points = get_ellipse_points(pos, V, B)

                    ! we have 2 points per ellipse
                    do p = 1, 2

                        ! ensure point is within the domain
                        call apply_periodic_bc(points(p, :))

                        ! get interpolation weights and mesh indices
                        call trilinear(points(p, :), is, js, weights)

                        do l = 1, ngp
                            volg(js(l), is(l)) = volg(js(l), is(l)) &
                                               + dble(m) * f12 * weights(l) * pvol
                        enddo
                    enddo
                enddo
            enddo
        end subroutine vol2grid_elliptic_symmetry_check


        subroutine vol2grid_non_elliptic
            integer          :: n, l
            double precision :: pos(2)
            double precision :: pvol

            do n = 1, n_parcels

                pos = parcels%position(n, :)
                pvol = parcels%volume(n)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call trilinear(pos, is, js, weights)

                do l = 1, ngp
                    volg(js(l), is(l)) = volg(js(l), is(l)) + weights(l) * pvol
                enddo
            enddo

        end subroutine vol2grid_non_elliptic


        subroutine par2grid
            vortg = zero
            volg = zero
            nparg = zero
            tbuoyg = zero

            if (parcel%is_elliptic) then
                call par2grid_elliptic
            else
                call par2grid_non_elliptic
            endif

            ! apply free slip boundary condition
            volg(0,  :) = two * volg(0,  :)
            volg(nz, :) = two * volg(nz, :)

            ! free slip boundary condition is reflective with mirror
            ! axis at the physical domain
            volg(1,    :) = volg(1,    :) + volg(-1,   :)
            volg(nz-1, :) = volg(nz-1, :) + volg(nz+1, :)

            vortg(0,  :) = two * vortg(0,  :)
            vortg(nz, :) = two * vortg(nz, :)
            vortg(1,    :) = vortg(1,    :) + vortg(-1,   :)
            vortg(nz-1, :) = vortg(nz-1, :) + vortg(nz+1, :)

            tbuoyg(0,  :) = two * tbuoyg(0,  :)
            tbuoyg(nz, :) = two * tbuoyg(nz, :)
            tbuoyg(1,    :) = tbuoyg(1,    :) + tbuoyg(-1,   :)
            tbuoyg(nz-1, :) = tbuoyg(nz-1, :) + tbuoyg(nz+1, :)

            ! exclude halo cells to avoid division by zero
            vortg(0:nz, :) = vortg(0:nz, :) / volg(0:nz, :)
            tbuoyg(0:nz, :) = tbuoyg(0:nz, :) / volg(0:nz, :)

            ! sum halo contribution into internal cells
            ! (be aware that halo cell contribution at upper boundary
            ! are added to cell nz)
            nparg(0,    :) = nparg(0,    :) + nparg(-1, :)
            nparg(nz-1, :) = nparg(nz-1, :) + nparg(nz, :)

            ! sanity check
            if (sum(nparg(0:nz-1, :)) /= n_parcels) then
                print *, "par2grid: Wrong total number of parcels!"
                stop
            endif

        end subroutine par2grid

        subroutine par2grid_elliptic
            double precision  :: points(2, 2)
            integer           :: n, p, l, i, j
            double precision  :: pvol, pvor, weight

            do n = 1, n_parcels
                pvol = parcels%volume(n)

                points = get_ellipse_points(parcels%position(n, :), &
                                            pvol, parcels%B(n, :))

                call get_index(parcels%position(n, :), i, j)
                i = mod(i + nx, nx)
                nparg(j, i) = nparg(j, i) + 1

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(p, :), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    ! the weight is halved due to 2 points per ellipse
                    do l = 1, ngp

                        weight = f12 * weights(l) * pvol

                        vortg(js(l), is(l)) = vortg(js(l), is(l)) &
                                            + weight * parcels%vorticity(n)

                        tbuoyg(js(l), is(l)) = tbuoyg(js(l), is(l)) &
                                             + weight * parcels%buoyancy(n)

                        volg(js(l), is(l)) = volg(js(l), is(l)) &
                                           + weight
                    enddo
                enddo
            enddo
        end subroutine par2grid_elliptic


        subroutine par2grid_non_elliptic
            integer          :: n, l, i, j
            double precision :: pos(2)
            double precision :: pvol, weight

            do n = 1, n_parcels

                pos = parcels%position(n, :)
                pvol = parcels%volume(n)

                call get_index(pos, i, j)
                i = mod(i + nx, nx)
                nparg(j, i) = nparg(j, i) + 1

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call trilinear(pos, is, js, weights)

                ! loop over grid points which are part of the interpolation
                do l = 1, ngp

                    weight = pvol * weights(l)

                    ! the weight is halved due to 2 points per ellipse
                    vortg(js(l), is(l)) = vortg(js(l), is(l))  &
                                        + weight * parcels%vorticity(n)

                    tbuoyg(js(l), is(l)) = tbuoyg(js(l), is(l)) &
                                         + weight * parcels%buoyancy(n)

                    volg(js(l), is(l)) = volg(js(l), is(l)) + weight
                enddo
            enddo

        end subroutine par2grid_non_elliptic


        subroutine grid2par(vel, vor, vgrad)
            double precision,       intent(inout) :: vel(:, :), vor(:), vgrad(:, :)

            if (parcel%is_elliptic) then
                   call grid2par_elliptic(vel, vor, vgrad)
            else
                   call grid2par_non_elliptic(vel, vor, vgrad)
            endif

        end subroutine grid2par


        subroutine grid2par_add(vel, vor, vgrad)
            double precision,       intent(inout) :: vel(:, :), vor(:), vgrad(:, :)

            if (parcel%is_elliptic) then
                   call grid2par_elliptic(vel, vor, vgrad, add=.true.)
            else
                   call grid2par_non_elliptic(vel, vor, vgrad, add=.true.)
            endif

        end subroutine grid2par_add


        subroutine grid2par_elliptic(vel, vor, vgrad, add)
            double precision,     intent(inout) :: vel(:, :), vor(:), vgrad(:, :)
            logical, optional, intent(in)       :: add
            integer                             :: ncomp
            double precision                    :: points(2, 2), weight
            integer                             :: n, p, c, l

            ! number of field components
            ncomp = 2

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                    vel(1:n_parcels, :) = zero
                    vor(1:n_parcels)    = zero
               endif
            else
               vel(1:n_parcels, :) = zero
               vor(1:n_parcels)    = zero
            endif

            vgrad(1:n_parcels, :) = zero

            do n = 1, n_parcels

                points = get_ellipse_points(parcels%position(n, :), &
                                            parcels%volume(n),      &
                                            parcels%B(n, :))

                ! we have 2 points per ellipse
                do p = 1, 2

                    ! ensure point is within the domain
                    call apply_periodic_bc(points(p, :))

                    ! get interpolation weights and mesh indices
                    call trilinear(points(p, :), is, js, weights)

                    ! loop over grid points which are part of the interpolation
                    do l = 1, ngp
                        weight = f12 * weights(l)

                        ! loop over field components
                        do c = 1, ncomp
                            ! the weight is halved due to 2 points per ellipse
                            vel(n, c) = vel(n, c) &
                                      + weight * velog(js(l), is(l), c)
                        enddo

                        do c = 1, 4
                            vgrad(n, c) = vgrad(n, c) &
                                        + weight * velgradg(js(l), is(l), c)
                        enddo

                        vor(n) = vor(n) + weight * vtend(js(l), is(l))
                    enddo
                enddo
            enddo

        end subroutine grid2par_elliptic


        subroutine grid2par_non_elliptic(vel, vor, vgrad, add)
            double precision,     intent(inout) :: vel(:, :), vor(:), vgrad(:, :)
            logical, optional, intent(in)       :: add
            integer                             :: ncomp
            integer                             :: n, c, l
            double precision                    :: pos(2)

            ! number of field components
            ncomp = 2

            ! clear old data efficiently
            if(present(add)) then
               if(add .eqv. .false.) then
                   vel(1:n_parcels, :) = zero
                   vor(1:n_parcels)    = zero
               endif
            else
               vel(1:n_parcels, :) = zero
               vor(1:n_parcels)    = zero
            endif

            vgrad(1:n_parcels, :) = zero

            do n = 1, n_parcels

                pos = parcels%position(n, :)

                ! ensure parcel is within the domain
                call apply_periodic_bc(pos)

                ! get interpolation weights and mesh indices
                call trilinear(pos, is, js, weights)

                ! loop over grid points which are part of the interpolation
                do l = 1, ngp
                    ! loop over field components
                    do c = 1, ncomp
                        vel(n, c) = vel(n, c) &
                                  + weights(l) * velog(js(l), is(l), c)
                    enddo

                    do c = 1, 4
                        vgrad(n, c) = vgrad(n, c) &
                                    + weights(l) * velgradg(js(l), is(l), c)
                    enddo

                    vor(n) = vor(n) &
                           + weights(l) * vtend(js(l), is(l))

                enddo
            enddo

        end subroutine grid2par_non_elliptic


        !
        ! tri-linear interpolation
        !
        subroutine trilinear(pos, ii, jj, ww)
            double precision, intent(in)  :: pos(2)
            integer,          intent(out) :: ii(4), jj(4)
            double precision, intent(out) :: ww(4)
            double precision              :: xy(2)

            ! (i, j)
            call get_index(pos, ii(1), jj(1))
            call get_position(ii(1), jj(1), xy)
            ww(1) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j)
            ii(2) = ii(1) + 1
            jj(2) = jj(1)
            call get_position(ii(2), jj(2), xy)
            ww(2) = product(one - abs(pos - xy) * dxi)

            ! (i, j+1)
            ii(3) = ii(1)
            jj(3) = jj(1) + 1
            call get_position(ii(3), jj(3), xy)
            ww(3) = product(one - abs(pos - xy) * dxi)

            ! (i+1, j+1)
            ii(4) = ii(2)
            jj(4) = jj(3)
            call get_position(ii(4), jj(4), xy)
            ww(4) = product(one - abs(pos - xy) * dxi)

            ! account for x periodicity
            call periodic_index_shift(ii)

        end subroutine trilinear

end module parcel_interpl
