module field_bc
    use parameters, only : mesh
    implicit none

    contains

        subroutine apply_field_bc(field)
            double precision, intent(inout) :: field(0:, 0:, :)
            integer                         :: i

            ! sum halo cell values to internal cells at the boundary
            field(1, :, :)            = field(1, :, :) + field(0, :, :)
            field(mesh%grid(1), :, :) = field(mesh%grid(1), :, :)      &
                                      + field(mesh%grid(1) + 1, :, :)

            field(:, 1, :)            = field(:, 1, :) + field(:, 0, :)
            field(:, mesh%grid(2), :) = field(:, mesh%grid(2), :)      &
                                      + field(:, mesh%grid(2) + 1, :)

            ! vertical bc
            call apply_free_slip(field)

            ! horizontal bc
            call apply_periodic(field)

        end subroutine apply_field_bc


        subroutine apply_periodic(field)
            double precision, intent(inout) :: field(0:, 0:, :)

            field(1, :, :)            = field(1, :, :)              &
                                      + field(mesh%grid(1), :, :)
            field(mesh%grid(1), :, :) = field(1, :, :)
        end subroutine apply_periodic


        subroutine apply_free_slip(field)
            double precision, intent(inout) :: field(0:, 0:, :)

            field(:, 1, :)            = 2.0 * field(:, 1, :)
            field(:, mesh%grid(2), :) = 2.0 * field(:, mesh%grid(2), :)
        end subroutine apply_free_slip

end module field_bc
