! =============================================================================
!               Permutation module implemented according to
!               the non-recursive Heap algorithm
!               (https://en.wikipedia.org/wiki/Heap%27s_algorithm)
! =============================================================================
module permute
    implicit none

    private

    public :: permute_generate

    contains

        subroutine permute_alloc(permutes, n)
            integer, allocatable, intent(inout) :: permutes(:, :)
            integer, intent(in) :: n
            integer             :: i, p

            p = 1
            do i = 1, n
                p = p * i
            enddo

            allocate(permutes(p, n))

        end subroutine permute_alloc

        ! non-recursive Heap algorithm
        ! Implemented according to
        ! https://en.wikipedia.org/wiki/Heap%27s_algorithm
        ! (visited 10 August 2021)
        subroutine permute_generate(permutes, n)
            integer, allocatable, intent(inout) :: permutes(:, :)
            integer                             :: array(0:n-1)
            integer, intent(in)                 :: n
            integer                             :: tmp, i, k
            integer                             :: c(n)

            call permute_alloc(permutes, n)

            do i = 0, n-1
                array(i) = i + 1
            enddo

            c = 0

            k = 1
            permutes(1, :) = array

            i = 1
            do while (i < n)
                if (c(i) < i) then
                    if (mod(i, 2) == 0) then
                        tmp = array(0)
                        array(0) = array(i)
                        array(i) = tmp
                    else
                        tmp = array(c(i))
                        array(c(i)) = array(i)
                        array(i) = tmp
                    endif
                    k = k + 1
                    permutes(k, :) = array

                    c(i) = c(i) + 1
                    i = 1
                else
                    c(i) = 0
                    i = i + 1
                endif
            enddo
        end subroutine permute_generate
end module permute
