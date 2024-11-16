module legendre_polynoms
! This module is used to prepare array of roots of Legendre polynomials
! for Gaussian quadrature function
use :: precision
implicit none

contains

function polynomial_division (f, g, r_in) result (q)
! Divinding of polynomial f, such as deg(f) = n >= 0, by polynomial
! g, such as deg(g) = m >= 0
!
! Input: polynoms f and g in the shape of real arrays of 
! their coefficients; optional pointer to polynom r_in (remainder of the division)
!
! Output: polynom q, deg(q) = n-m, result of division; pointer to 
! polynom r_in

        real(prec), dimension (0:) :: f
        real(prec), dimension (0:) :: g
        real(prec), dimension (0:size(f)-1), optional :: r_in

        real(prec), dimension (0:size(f) - size(g)) :: q
        real(prec), dimension (0:size(f) - 1) :: r

        integer :: i, n, m
        real(prec), dimension (0:size(f) - 1) :: intermediate

        n = size(f) - 1
        m = size(g) - 1

        r = f

        ! Деление в столбик
        do i = 0, n-m
               intermediate = 0 
               intermediate(i : i+m) = g*r(i)

               q(i) = r(i)
               r = r - intermediate 
        enddo 


        if (present(r_in)) r_in = r
        
end function polynomial_division

function make_legendre (n) result (p)
! Makes Legendre polynomial of degree n
! using iterative algorithm
! 
! Input: integer degree of polynomial n
!
! Output: array of coefficients of legendre polynomial (real(prec) (0:n))

        integer :: n
        real(prec), dimension (0:n) :: p

        integer :: i, k
        real(prec) :: nr
        real(prec), dimension (0:n) :: p1, p2

        p1 = 0.0_prec
        p2 = 0.0_prec

        if (n > 1) then
                p2(n) = 1.0_prec
                p1(n-1) = 1.0_prec
        endif

        do i = 2, n
                do k = 0, n-1
                        nr = real(i, prec) 

                        p(k) = (2.0_prec*nr - 1.0_prec)*&
                                p1(k+1)/nr - (nr-1.0_prec) * p2(k)/nr
                enddo
                        p(n) = -(nr-1.0_prec) * p2(n)/nr

                p2 = p1
                p1 = p
        enddo 

        select case (n)
                case (0)
                        p = 1.0_prec
                case (1)
                        p(1) = 0.0_prec
                        p(0) = 1.0_prec
                case default
                        p = p1
        end select

end function make_legendre

function horner (x, pol) result (res)
! Calculates polynomial value at point x
! 
! Input: point x (real(prec)); polynomial pol (real(prec) (0:n)) as array of its contains  
! (pol[0]*x^n + ... + pol[n])
!
! Output: value of polynomial at point x res (real(prec))

        real(prec) :: x, res
        real(prec), dimension(0:) :: pol

        integer :: n, i

        n = size(pol)-1

        res = 0
        do i = 0, n-1
                res = (res + pol(i)) * x
        enddo
                res = res + pol(n)
end function horner

function bernulli_method(polynom) result (roots)
! Finds roots of legendre polynomial using Bernulli method
! 
! Input: legendre polynomial as array of its coefficients (real(prec) (0:n))
! 
! Output: array of roots of legendre polynomial (real (1:n))

        real(prec), dimension (0:) :: polynom
        real(prec), dimension (1:size(polynom)-1) :: roots

        integer :: n, i

        real(prec), dimension (:), allocatable :: y
        real(prec), dimension (0:size(polynom)-1) :: temp

        real(prec) :: yn
        real(prec), dimension (1:2) :: x_history

        n = size(polynom) - 1

        ! Если полином первой или второй степени, то ответ выражается явно
        ! (алгоритм, описанный далее не работает)
        if (n == 1) then
                roots = -polynom(1)/polynom(0)
                return
        endif
        if (n == 2) then
                roots(1) = sqrt(-polynom(2)/polynom(0))
                roots(2) = -sqrt(-polynom(2)/polynom(0))
                return
        endif


        ! Присвоение начального значения полинома, чьи корни вычисляются
        ! (деление производится для уменьшения количества действий при 
        ! реализации метода Бернулли)
        temp = -polynom/polynom(0)

        ! Поиск каждого второго корня полинома (остальные симметричны)
        ! 
        ! Если степень полинома четная, то ищутся все корни, кроме двух,
        ! в обратном случае - все, кроме одного
        do i = 1, n - mod(n+1, 2) - 1, 2
                allocate (y(0:n-i))
                call random_number (y)

                ! Первое приближение корня
                x_history(2) = 0
                x_history(1) = y(0)/y(1)

                ! Поиск корня производится до момента, когда следующий 
                ! корень мало отличается от предыдущего
                do while (abs(x_history(1) - x_history(2)) > ZERO)

                        ! Два шага метода Бернулли, т.к. последовательность не сходится
                        ! ввиду симметричности корней полинома Лежандра
                        yn = dot_product(temp(i:n), y)
                                y(1:n-i) = y(0:n-i-1) 
                                y(0) = yn
                        yn = dot_product(temp(i:n), y)
                                y(1:n-i) = y(0:n-i-1) 
                                y(0) = yn

                        ! Нахождение нового приблтжения корня (квадрата корня) 
                        x_history(2) = x_history(1)
                        x_history(1) = y(0)/y(2)

                enddo


                ! Деление полинома. Деление производится с целью исключения найденных корней из него
                temp(i+1:n) = polynomial_division(temp(i-1:n), [1.0_prec, 0.0_prec, -abs(x_history(1))])
                temp(i-1) = 0
                temp(i) = 0


                ! Заолнение массива корней
                roots(i) = sqrt(abs(x_history(1)))
                roots(i+1) = -sqrt(abs(x_history(1)))

                deallocate(y)
        enddo
        
        ! Явное выражение последних корней полинома
        ! (для них алгоритм не раюботает)
        if (mod(n, 2) == 1) roots(n) = -temp(n)/temp(n-1)
        if (mod(n, 2) == 0) then
                roots(n-1) = sqrt(abs(temp(n)/temp(n-2)))
                roots(n) = -sqrt(abs(temp(n)/temp(n-2)))
        endif

end function bernulli_method

function roots_legendre(n) result (roots)
! Returns roots of legendre polynomial of degree n
! 
! Input: integer degree of legendre polynomial n
!
! Output: array of roots of legendre polynomial (real (1:n))

        integer :: n
        real(prec), dimension (1:n) :: roots

        real(prec), dimension (0:n) :: legangre


        legangre = make_legendre(n)
        roots = bernulli_method(legangre)
end function roots_legendre

end module legendre_polynoms
