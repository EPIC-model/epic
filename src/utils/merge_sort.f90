! =============================================================================
!                                   Mergesort
!   Sort an array of values in ascending order and return also the
!   re-ordering of the indices (in order to sort other arrays accordingly).
!
!   Implemented according to Ottmann, T. and Widmayer, P., Algorithmen
!   und Datenstrukturen, 4. Auflage, 2002, Spektrum Akademischer Verlag
! =============================================================================
module merge_sort
    implicit none

    private :: mergesort, merging, imergesort, imerging

    ! 23 August 2021
    ! https://fortran-lang.discourse.group/t/question-about-interfaces-and-submodules/584
    interface msort
        module procedure :: dmsort
        module procedure :: imsort
    end interface

    contains

        ! a  : value array to sort
        ! ai : index array
        subroutine dmsort(a, ai)
            double precision, intent(inout) :: a(:)
            integer, intent(out)            :: ai(:)
            integer                         :: n, m

            m = size(a)

            do n = 1, m
                ai(n) = n
            enddo

            call mergesort(a, ai, 1, m)

        end subroutine dmsort

        recursive subroutine mergesort(a, ai, l, r)
            double precision, intent(inout) :: a(:)
            integer,          intent(inout) :: ai(:)
            integer,          intent(in)    :: l, r
            integer                         :: m

            if (l < r) then
                m = (l + r) / 2

                call mergesort(a, ai, l, m)
                call mergesort(a, ai, m + 1, r)

                call merging(a, ai, l, m, r)
            endif
        end subroutine mergesort

        subroutine merging(a, ai, l, m, r)
            double precision, intent(inout) :: a(:)
            integer,          intent(inout) :: ai(:)
            integer,          intent(in)    :: l, m, r
            double precision                :: b(r)
            integer                         :: bi(r)
            integer                         :: i, j, k, h

            i = l
            j = m + 1
            k = l

            do while ((i <= m) .and. (j <= r))
                if (a(i) <= a(j)) then
                    b(k) = a(i)
                    bi(k) = ai(i)
                    i = i + 1
                else
                    b(k) = a(j)
                    bi(k) = ai(j)
                    j = j + 1
                endif
                k = k + 1
            enddo

            if (i > m) then
                do h = j, r
                    b(k + h - j) = a(h)
                    bi(k + h - j) = ai(h)
                enddo
            else
                do h = i, m
                    b(k + h - i) = a(h)
                    bi(k + h - i) = ai(h)
                enddo
            endif

            do h = l, r
                a(h) = b(h)
                ai(h) = bi(h)
            enddo
        end subroutine merging



        ! a  : value array to sort
        ! ai : index array
        subroutine imsort(a, ai)
            integer, intent(inout) :: a(:)
            integer, intent(out)            :: ai(:)
            integer                         :: n, m

            m = size(a)

            do n = 1, m
                ai(n) = n
            enddo

            call imergesort(a, ai, 1, m)

        end subroutine imsort

        recursive subroutine imergesort(a, ai, l, r)
            integer, intent(inout) :: a(:)
            integer, intent(inout) :: ai(:)
            integer, intent(in)    :: l, r
            integer                :: m

            if (l < r) then
                m = (l + r) / 2

                call imergesort(a, ai, l, m)
                call imergesort(a, ai, m + 1, r)

                call imerging(a, ai, l, m, r)
            endif
        end subroutine imergesort

        subroutine imerging(a, ai, l, m, r)
            integer, intent(inout) :: a(:)
            integer, intent(inout) :: ai(:)
            integer, intent(in)    :: l, m, r
            integer                :: b(r)
            integer                :: bi(r)
            integer                :: i, j, k, h

            i = l
            j = m + 1
            k = l

            do while ((i <= m) .and. (j <= r))
                if (a(i) <= a(j)) then
                    b(k) = a(i)
                    bi(k) = ai(i)
                    i = i + 1
                else
                    b(k) = a(j)
                    bi(k) = ai(j)
                    j = j + 1
                endif
                k = k + 1
            enddo

            if (i > m) then
                do h = j, r
                    b(k + h - j) = a(h)
                    bi(k + h - j) = ai(h)
                enddo
            else
                do h = i, m
                    b(k + h - i) = a(h)
                    bi(k + h - i) = ai(h)
                enddo
            endif

            do h = l, r
                a(h) = b(h)
                ai(h) = bi(h)
            enddo
        end subroutine imerging

end module merge_sort
