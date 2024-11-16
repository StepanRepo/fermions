module gauss_int
        use :: sle_solve
        use :: legendre_polynoms
        use :: precision
implicit none

contains

subroutine get_coefficients (t) 
! Calculates coefficients of gaussian quadrature
!
! Input: 
! t -- real(prec)(1:n) -- array of roots of legendre polynomials 
!
! Output: file in directory "./quads/" entitled "quad??.dat", where ?? is n
! File format: coefficient[1] t[1]
!                        ...
!               coefficient[n] t[n]

        integer, parameter :: output_file = 2
        character(18) :: file_name
        logical :: directory_exists

        real(prec), dimension(1:) :: t

        real(prec), dimension(1:size(t), 1:size(t)) :: matrix
        real(prec), dimension(1:size(t)) :: coefficients, B 

        integer :: i, n

        n = size(t)

        ! Составление матрицы системыт уравнений для поиска коэффициентов
        ! Вычисление вектора свободных членов
        do i = 1, n
                matrix(i, :) = t(:)**real(i-1, prec)
                B(i) = 2.0_prec*real(mod(i, 2), prec)/real(i, prec)
        enddo

        ! Нахождение коэффициентов с помощью решения СЛУ
        ! методом Гаусса с выбором главного элемента
        coefficients = choice(matrix, B)

        ! Составление имени файла
        write(file_name, "('./quads/quad', i2.2, '.dat')") n

        ! Check if directory exists
        ! make one if not
        inquire(file = "quads", exist = directory_exists)

        if (.not. directory_exists) &
                call execute_command_line("mkdir quads")

        ! Вывод найденных коэффициентов в файл
        open (output_file, file = file_name)
        do i = 1, n
                write(output_file, *) coefficients(i), t(i)
        enddo

        close(output_file)
end subroutine get_coefficients



function gint (f, a, b, n_in) result (res)
! Calculates integral of function f with limits a and b using Gaussian quadrature
! (coefficients for quadrature are readed from file "./quads/quad??.dat" where ?? is n)
! 
! Input: 
!       f -- name of function, which gets real(prec) input value and returns real(prec) output value;

!       a -- real(prec) -- lower limit of integration  
!       b -- real(prec) -- upper limit of integration 
!       n -- integer -- optional, degree of fittet polynomial
! 
! Output: 
!       res -- real(prec) -- estimated value of the integral 

        interface
                function f(x) result (res)
                        use :: precision
                        real(prec) :: x, res
                end function f
        end interface

        real(prec) :: a, b
        integer, optional :: n_in
        real(prec) :: res

        real(prec) :: t
        real(prec), dimension (:), allocatable :: coefficients, ft
        integer :: n, i

        character(18) :: file_name
        logical :: file_exists

        ! Стандартное количество разбиений - 5
        n = 5
        if (present(n_in)) n = n_in


        ! Нахождение имени файла
        write(file_name, "('./quads/quad', i2.2, '.dat')") n

        ! Check if file exists. If not, calculate the coefficients
        inquire(file = file_name, exist = file_exists)
        
        if (.not. file_exists) &
                call get_coefficients(roots_legendre(n))


        ! Чтение коэффициентов квадратуры и корней полинома Лежандра
        open(1, file = file_name)
        allocate (coefficients(1:n), ft(1:n))
              
        do i = 1, n
                read (1, *) coefficients(i), t

                ! Вычисление значений функции в точках, связанных с корнями
                ft(i) = f((b + a + (b - a)*t)/2.0_prec)
        enddo

        ! Интеграл от функции
        res = dot_product(coefficients, ft)*(b - a)/2.0_prec

        close (1)


end function gint

function gint_to_inf (f, a, step, n_in) result (res)
        interface
                function f(x) result (res)
                        use :: precision
                        real(prec) :: x, res
                end function f
        end interface

        real(prec) :: a, step
        integer, optional :: n_in
        real(prec) :: res

        integer :: n, i
        real(prec) :: x0, x1, partial

        ! Стандартное количество разбиений - 5
        n = 5
        if (present(n_in)) n = n_in
        

        res = 0.0_prec
        i = 0

        x0 = 0.0_prec
        x1 = step


        do while (abs(partial) > ZERO .or. i < 2)

                ! calculate one iteration of the integral
                partial = gint (f, x0, x1, 5)
                res = res + partial

                ! reset iteraion vars to new cycle
                x0 = x1
                x1 = x1 + step

                i = i+1

                if (i > 1000) exit
        enddo


        if (i < 3) then
                write(*, *) "Integral conberges too fast"
        endif

        if (i > 1000) then
                write(*, *) "Integral is not convergent after 10^3 iterations"
        endif

end function gint_to_inf

        


end module gauss_int
