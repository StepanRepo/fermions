#include "libraries.f95"

module sle_solve
        use :: precision
implicit none

interface solve
        module procedure :: choice
end interface solve

public :: solve
private :: LAPACK_wrapper

contains

function LAPACK_wrapper (A_in, B_in) result (X)
        real(prec) :: A_in(1:, 1:), B_in(1:)
        real(prec) :: X(1:size(B_in))

        real(prec) :: A(1:size(B_in), 1:size(B_in)), B(1:size(B_in), 1)
        integer :: pivot(size(B_in)), stat, n

        n = size(B_in)

        A = A_in
        B(:, 1) = B_in

        if (prec == 4) then
                call sgesv(n, 1, A, n, pivot, B, n, stat)
        else if (prec == 8) then
                call dgesv(n, 1, A, n, pivot, B, n, stat)
        endif

        if (stat > 0) then
                write(*, *) "Can't solve"
        endif

        X = B(:, 1)
end function

function choice (A, B) result (X)
        real(prec) :: A(1:, 1:), B(1:)
        real(prec) :: X(1:size(B)), temp(1:size(B))

        integer :: i, j, n

        integer :: ml(2), transforms(2, size(B))
        ! tansforms:
        ! | old row | old col | new row | new col |
        ! |   1     |   1     |   ...   |   ...   |
        ! |   2     |   2 !    |   ...   |   ...   |
                          !
        n = size (B)      !
                          !

        if (IS_LAPACK .and. prec <= 8) then
                X = LAPACK_wrapper(A, B)
                return
        endif

        do i = 1, n

                ! Выбор максимального элемента
                ml = maxloc(abs(A(i:n, i:n))) + i - 1 ! Добавка из-за сокращения матрицы на входе

                temp = A(i, :)
                        A(i, :) = A(ml(1), :)
                        A(ml(1), :) = temp

                temp(1) = B(i)
                        B(i) = B(ml(1))
                        B(ml(1)) = temp(1)

                temp = A(:, i)
                        A(:, i) = A(:, ml(2))
                        A(:, ml(2)) = temp

                transforms(:, i) = ml   ! Запомнить перестановку


                if (abs(A(i,i)) < ZERO) write(*, "('Деление на ноль!')") 
                        ! Проверка деления на близкон к нулю число


                B(i) = B(i) / A(i, i)           ! Деление свободного вектора до переопределения знаменателя
                                                ! иначе это просто деление на 1
                A(i, :) = A(i, :) / A(i, i)

                do j = 1, n
                        if (j == i) cycle
                        
                        B(j) = B(j) - B(i)*A(j, i)
                        A(j, :) = A(j, :) - A(i, :)*A(j, i)
                enddo 
        enddo

        ! Расстановка всех элементов на прежние места
        do i = n, 1, -1
                ml = transforms(:, i)

                temp = A(i, :)
                        A(i, :) = A(ml(1), :)
                        A(ml(1), :) = temp

                temp(1) = B(i)
                        B(i) = B(ml(1))
                        B(ml(1)) = temp(1)

                temp = A(:, i)
                        A(:, i) = A(:, ml(2))
                        A(:, ml(2)) = temp
        enddo

        ! Заполнение тех элементов Х, чей коэффициент в строке не равен 0 (в данном случае масимален)
        do i = 1, n
                X(maxloc(A(i,:))) = B(i)
        enddo

end function choice

end module sle_solve
