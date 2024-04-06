! =============================================================================
!       Parcel ellipse container for the 2D surface parcels.
!       Note that at the surfaces the flow is compressible.
! =============================================================================
module parcel_ellipse
    use parcel_container, only : pc_type
    use constants, only : pi   &
                        , fpi  &
                        , zero &
                        , two  &
                        , f12  &
                        , f14
    use mpi_utils, only : mpi_exit_on_error
    use parameters, only : amin
    use armanip, only : resize_array
    implicit none

    type, extends(pc_type) :: ellipse_pc_type

        double precision :: z_position = 0.0d0 !FIXME

        ! ---------------------------------------------------------------------
        ! Additional parcel attributes:
        double precision, allocatable, dimension(:) :: area
        double precision, allocatable, dimension(:) :: circ  ! circulation

        integer :: IDX_AREA
        integer :: IDX_CIRC

        contains

            procedure          :: allocate => parcel_ellipse_allocate
            procedure          :: deallocate => parcel_ellipse_deallocate
            procedure          :: replace => parcel_ellipse_replace
            procedure          :: resize => parcel_ellipse_resize
            procedure          :: serialize => parcel_ellipse_serialize
            procedure          :: deserialize => parcel_ellipse_deserialize
            procedure          :: get_points => parcel_ellipse_get_points
            procedure          :: get_eigenvalue => parcel_ellipse_get_eigenvalue
            procedure          :: get_eigenvector => parcel_ellipse_get_eigenvector
            procedure          :: get_angle => parcel_ellipse_get_angle
            procedure          :: get_ab => parcel_ellipse_get_ab
            procedure          :: get_area => parcel_ellipse_get_area
            procedure          :: get_aspect_ratio => parcel_ellipse_get_aspect_ratio
            procedure          :: is_small => parcel_ellipse_is_small
            procedure          :: get_z_position => parcel_ellipse_get_z_position
            procedure          :: get_position => parcel_ellipse_get_position
            procedure          :: set_position => parcel_ellipse_set_position
            procedure          :: set_volume => parcel_ellipse_set_area
            procedure          :: get_dBdt => parcel_ellipse_get_dBdt


            !------------------------------------------------------------------
            ! Procedures in submodules:
            procedure :: init => parcel_ellipse_init
            procedure :: split => parcel_ellipse_split

    end type ellipse_pc_type

    !--------------------------------------------------------------------------
    ! Define interface for submodule routines:
    interface
        ! Implemented in parcel_ellipse_init
        module subroutine parcel_ellipse_init(this)
            class(ellipse_pc_type), intent(inout) :: this
        end subroutine parcel_ellipse_init

        module subroutine parcel_ellipse_split(this)
            class(ellipse_pc_type), intent(inout) :: this
        end subroutine parcel_ellipse_split

    end  interface

    contains

        ! Sets the buffer indices for sending parcels around
        subroutine parcel_ellipse_allocate(this, num)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: num
            integer                               :: i

            call this%parcel_base_allocate(num=num,    &
                                           n_pos=2,    &
                                           n_vec=3,    &
                                           n_shape=3,  &   ! the flow is compressible
                                           n_strain=4)

            allocate(this%area(num))
            allocate(this%circ(num))

            this%IDX_POS_BEG   = 1   ! x-position
            this%IDX_POS_END   = 2   ! y-position
            this%IDX_VOR_BEG   = 3   ! x-vorticity
            this%IDX_VOR_END   = 4   ! z-vorticity
            this%IDX_SHAPE_BEG = 5   ! B11 shape matrix element
            this%IDX_SHAPE_END = 7   ! B22 shape matrix element
            this%IDX_AREA      = 8   ! area
            this%IDX_CIRC      = 9
            this%IDX_VOL       = 10  ! volume (needed for mixing with interior parcels)
            this%IDX_BUO       = 11  ! buoyancy

            i = this%IDX_BUO + 1
#ifndef ENABLE_DRY_MODE
            this%IDX_HUM  = i
            i = i + 1
#endif

            ! LS-RK variables
            this%IDX_RK_POS_BEG    = i
            this%IDX_RK_POS_END    = i + 1
            this%IDX_RK_VOR_BEG    = i + 2
            this%IDX_RK_VOR_END    = i + 3
            this%IDX_RK_SHAPE_BEG  = i + 4
            this%IDX_RK_SHAPE_END  = i + 6
            this%IDX_RK_STRAIN_BEG = i + 7
            this%IDX_RK_STRAIN_END = i + 10

            this%attr_num = this%IDX_RK_STRAIN_END

            ! these parcels only are in one vertical grid cell
            this%nz = 1

        end subroutine parcel_ellipse_allocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_deallocate(this)
            class(ellipse_pc_type), intent(inout) :: this

            call this%parcel_base_deallocate

            if (.not. allocated(this%area)) then
                return
            endif

            deallocate(this%area)
            deallocate(this%circ)

        end subroutine parcel_ellipse_deallocate

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_replace(this, n, m)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: n, m

            ! Call parent class subroutine
            call this%parcel_base_replace(n, m)

            this%area(n) = this%area(m)
            this%circ(n) = this%circ(m)

        end subroutine parcel_ellipse_replace

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_resize(this, new_size)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: new_size

            ! Call parent class subroutine
            call this%parcel_base_resize(new_size)

            call resize_array(this%area, new_size, n_copy=this%local_num)
            call resize_array(this%circ, new_size, n_copy=this%local_num)

        end subroutine parcel_ellipse_resize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_serialize(this, n, buffer)
            class(ellipse_pc_type), intent(in)  :: this
            integer,                intent(in)  :: n
            double precision,       intent(out) :: buffer(this%attr_num)

            ! Call parent class subroutine
            call this%parcel_base_serialize(n, buffer)

            buffer(this%IDX_AREA) = this%area(n)
            buffer(this%IDX_CIRC) = this%circ(n)

        end subroutine parcel_ellipse_serialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_deserialize(this, n, buffer)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: n
            double precision,       intent(in)    :: buffer(this%attr_num)

            ! Call parent class subroutine
            call this%parcel_base_deserialize(n, buffer)

            this%area(n) = buffer(this%IDX_AREA)
            this%circ(n) = buffer(this%IDX_CIRC)

        end subroutine parcel_ellipse_deserialize

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the ellipse support points for par2grid and grid2par
        ! @param[in] n parcel index
        ! @returns the parcel support points
        function parcel_ellipse_get_points(this, n) result(points)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: c, a2, evec(2), h(2)
            double precision                   :: points(2, 2)

            a2 = this%get_eigenvalue(n)

            c = dsqrt(dabs(two * a2 - this%B(1, n) - this%B(3, n)))

            evec = this%get_eigenvector(n, a2)

            h = f12 * c * evec

            points(:, 1) = this%position(:, n) - h
            points(:, 2) = this%position(:, n) + h

        end function parcel_ellipse_get_points

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the largest eigenvalue (i.e. the semi-major axis squared [a**2])
        ! @param[in] n parcel index
        ! @returns the largest eigenvalue
        function parcel_ellipse_get_eigenvalue(this, n) result(a2)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: a2

            a2 = f12 * (this%B(1, n) + this%B(3, n)) &
               + dsqrt(f14 * (this%B(1, n) - this%B(3, n)) ** 2 + this%B(2, n) ** 2)

        end function parcel_ellipse_get_eigenvalue

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the eigenvector of the largest eigenvalue
        ! @param[in] n parcel index
        ! @param[in] a2 is the largest eigenvalue
        ! @returns the eigenvector
        function parcel_ellipse_get_eigenvector(this, n, a2) result(evec)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision,       intent(in) :: a2
            double precision                   :: evec(2)

            evec(1) = a2 - this%B(3, n)
            evec(2) = this%B(2, n)

            if (dabs(evec(1)) + dabs(evec(2)) == zero) then
                if (this%B(1, n) > this%B(3, n)) then
                    evec(1) = evec(1) + epsilon(evec(1))
                else
                    evec(2) = evec(2) + epsilon(evec(2))
                endif
            endif

            evec = evec / norm2(evec)

        end function parcel_ellipse_get_eigenvector

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! used in unit tests only
        ! @param[in] n parcel index
        function parcel_ellipse_get_angle(this, n) result(angle)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: a2, evec(2)
            double precision                   :: angle

            a2 = this%get_eigenvalue(n)

            evec = this%get_eigenvector(n, a2)

            angle = datan2(evec(2), evec(1))

        end function parcel_ellipse_get_angle

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the product of the semi-minor and semi-major axis.
        ! @param[in] area of the parcel
        ! @returns ab = area / pi
        elemental function parcel_ellipse_get_ab(this, area) result(ab)
            class(ellipse_pc_type), intent(in) :: this
            double precision,       intent(in) :: area
            double precision                   :: ab

            ab = area * fpi

        end function parcel_ellipse_get_ab

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! @param[in] n parcel index
        ! @returns the parcel area
        function parcel_ellipse_get_area(this, n) result(area)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: area

#ifndef NDEBUG
            if (this%B(1, n) * this%B(3, n) < this%B(2, n) ** 2) then
                call mpi_exit_on_error(&
                    "in parcel_ellipse::get_area: Determinant is negative. Unable to calculate parcel area.")
            endif
#endif
            area = pi * dsqrt(this%B(1, n) * this%B(3, n) - this%B(2, n) ** 2)

        end function parcel_ellipse_get_area

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Obtain the aspect ratio a/b of the parcel(s).
        ! @param[in] n parcel index
        ! @param[in] a2 is the largest eigenvalue
        ! @returns lam = a/b
        elemental function parcel_ellipse_get_aspect_ratio(this, n, a2) result(lam)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision,       intent(in) :: a2
            double precision                   :: lam

            lam = (a2 / this%area(n)) * pi

        end function parcel_ellipse_get_aspect_ratio

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        pure function parcel_ellipse_is_small(this, n) result(l_small)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            logical                            :: l_small

            l_small = (this%area(n) < amin)

        end function parcel_ellipse_is_small

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        pure function parcel_ellipse_get_z_position(this, n) result(z_pos)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: z_pos

            z_pos = this%z_position

        end function parcel_ellipse_get_z_position

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        pure function parcel_ellipse_get_position(this, n) result(pos)
            class(ellipse_pc_type), intent(in) :: this
            integer,                intent(in) :: n
            double precision                   :: pos(3)

            pos(1) = this%position(1, n)
            pos(2) = this%position(2, n)
            pos(3) = this%z_position

        end function parcel_ellipse_get_position

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_set_position(this, n, pos)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: n
            double precision,       intent(in)    :: pos(3)

            this%position(1, n) = pos(1)
            this%position(2, n) = pos(2)

        end subroutine parcel_ellipse_set_position

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        subroutine parcel_ellipse_set_area(this, n, vol)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                intent(in)    :: n
            double precision,       intent(in)    :: vol

            this%area(n) = vol

        end subroutine parcel_ellipse_set_area

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! Advance the B matrix.
        ! @param[in] n parcel index
        ! @returns the updated B matrix components (B11, B12 and B22) in Bout
        function parcel_ellipse_get_dBdt(this, n) result(Bout)
            class(ellipse_pc_type), intent(inout) :: this
            integer,                  intent(in)  :: n
            double precision                      :: Bout(3)

            ! B11 = 2 * (dudx * B11 + dudy * B12)
            Bout(1) = two * (this%strain(1, n) * this%B(1, n) + &
                             this%strain(2, n) * this%B(2, n))

            ! B12 = (du/dx + dv/dy) * B12 + dvdx * B11 + dudy * B22
            Bout(2) = (this%strain(1, n) + this%strain(4, n)) * this%B(2, n) + &
                       this%strain(3, n)                      * this%B(1, n) + &
                       this%strain(2, n)                      * this%B(3, n)

            ! B22 = 2 * (dvdx * B12 + dvdy * B22)
            Bout(3) = two * (this%strain(3, n) * this%B(2, n) + &
                             this%strain(4, n) * this%B(3, n))

        end function parcel_ellipse_get_dBdt

end module parcel_ellipse
