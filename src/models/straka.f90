module straka
    use physics
    use options, only : parcel_info
    use constants
    use parcel_container, only : parcels, n_parcels
    use parcel_split, only : split_ellipses
    use parameters, only : dx, vcell
    use ellipse, only : get_ab, get_B22, get_eigenvalue
    use fields, only : get_position
    implicit none

    contains

        subroutine straka_init
            integer          :: n
            double precision :: xc, xr, zc, zr, L, dT, dtheta, theta_0, lam

            ! in metres
            xc = zero
            xr = four
            zc = three
            zr = two

            ! reference potential temperature
            theta_0 = 300.0d0 ! Kelvin

            ! aspect ratio: lam = a / b
            lam = max(dx(2) / dx(1), dx(1) / dx(2))

            call straka_refine(lam)

            dT = zero
            do n = 1, n_parcels
                L = ((parcels%position(n, 1) - xc) / xr) ** 2 &
                  + ((parcels%position(n, 2) - zc) / zr) ** 2

                L = dsqrt(L)

                ! temperature perturbation
                if (L > one) then
                    dT = zero
                else
                    dT = -7.5d0 * (dcos(pi * L) + one)
                endif

                ! potential temperatur
                dtheta = dT * pi

                ! MPIC paper:
                ! liquid-water buoyancy is defined by b = g * (theta − theta_0) / theta_0
                ! (dtheta = theta - theta_0)
                parcels%buoyancy(n) = gravity * dtheta / theta_0
            enddo
        end subroutine straka_init

        subroutine straka_refine(lam)
            double precision, intent(inout) :: lam
            integer                         :: n
            double precision                :: B22, a2

            if (.not. parcel_info%is_elliptic) then
                return
            endif

            ! we need to initialize all parcel B11 first
            do n = 1, n_parcels
                parcels%B(n, 1) = lam * get_ab(parcels%volume(n))
            enddo

            ! do refining by splitting
            do while (lam >= parcel_info%lambda)
                call split_ellipses(parcels, parcel_info%lambda, parcel_info%vmaxfraction)

                B22 = get_B22(parcels%B(1, 1), zero, parcels%volume(1))
                a2 = get_eigenvalue(parcels%B(1, 1), zero, B22)
                lam = a2 / get_ab(parcels%volume(1))
            end do
        end subroutine straka_refine

end module
