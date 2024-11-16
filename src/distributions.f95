module distributions
use :: gauss_int
use :: config
use :: precision
implicit none

!-----------------------------
! Fundamental constants (source: Wikipedia)
!-----------------------------

real(prec), parameter :: c = 29979245800.0_prec         ! speed of light, cm/s
real(prec), parameter :: h = 6.62607015e-27_prec        ! Planck's constant, erg*s
real(prec), parameter :: kB = 1.380649e-16_prec         ! Boltzmann constant, erg/K

real(prec), parameter :: pi = 3.14159265359_prec        ! pi




contains


!-----------------------------
! Auxillary functions
!-----------------------------
        function Fermi (n_in, xi) result (res)
                implicit none

                integer :: n_in
                real(prec) :: xi, res

                real(prec) :: n
                real(prec) :: partial, step

                ! set initial values
                n = real(n_in, prec)
                
                if (xi < 1.0_prec) then
                        step = 1.0_prec
                else
                        step = xi/10.0_prec
                endif

                
                res = gint_to_inf(f, 0.0_prec, step, 5)

                contains

                function f (x)
                        real(prec) :: x, f

                        f = (x**n) / (exp(x - xi) + 1.0_prec)
                end function f

        end function Fermi


!-----------------------------
! Convertation functions
!-----------------------------

        ! unit to unitless
        function unit2ul (u) result (ul)
        ! transforms (T, mu) to (y, phi)
        ! u is ordered as (T_i, mu_i)
        ! ul is oredered as (y_i, phi_i)
        
                real(prec), dimension(0:, 0:) :: u
                real(prec) :: energy
                real(prec), dimension(0:(size(u(:, 0)) - 1), 0:1) :: ul

                if (m > 0.0_prec) then
                        energy = m*c*c
                else if (m < 0.0_prec) then
                        write(*, *) "Particle mass can't be negative, m = ", m
                else
                        write(*, *) "Mass is zero using frequency"
                        energy = h*nu
                endif


                ul(:, 0) =  energy / kB / u(:, 0)
                ul(:, 1) =  (u(:, 1) + m*c*c) / kB / u(:, 0)
        end function unit2ul


        ! unitless to unit
        function ul2unit (ul) result (u)
        ! transforms (y, phi) to (T, mu)
        ! ul is oredered as (y_i, phi_i)
        ! u is ordered as (T_i, mu_i)
                real(prec), dimension(0:, 0:) :: ul
                real(prec) :: energy
                real(prec), dimension(0: (size(ul(:, 0)) - 1), 0:1) :: u

                if (m > 0.0_prec) then
                        energy = m*c*c
                else if (m < 0.0_prec) then
                        write(*, *) "Particle mass can't be negative, m = ", m
                else
                        write(*, *) "Mass is zero using frequency"
                        energy = h*nu
                endif

                u(:, 0) =  energy / kB / ul(:, 0)
                u(:, 1) = kB * u(:, 0) * ul(:, 1) - energy
        end function ul2unit


        function distr2phys (distr) result (phys)
                real(prec), dimension(0:, 0:) :: distr
                real(prec), dimension(0: (size(distr(:, 0)) - 1), 0:2) :: phys

                real(prec) :: lc, energy, g
                real(prec) :: coef

                if (m > 0.0_prec) then
                        energy = m*c*c
                        lc = h/m/c
                else if (m < 0.0_prec) then
                        write(*, *) "Particle mass can't be negative, m = ", m
                else
                        write(*, *) "Mass is zero using frequency"
                        energy = h*nu
                        lc = c/nu
                endif

                g = 2.0_prec*s + 1.0_prec

                coef = 4.0_prec * pi * g / (lc**3.0_prec)

                phys(:, 0) = coef * distr(:, 0)
                phys(:, 1) = coef * energy * distr(:, 1)
                phys(:, 2) = coef * energy / 3.0_prec * distr(:, 2)

        end function distr2phys



!-----------------------------
! Distribution functions for high temperature regime
! y -> 0+
!-----------------------------
        
        function F_n_high_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi

                xi = phi - y

                res = Fermi(2, xi) + 2.0_prec * y * Fermi(1, xi) + .5_prec*y*y * Fermi(0, xi)
                res = res / (y**3)

                res = res + exp(xi) / (exp(xi) + 1.0_prec) / 6.0_prec
                res = res - y/6.0_prec/16.0_prec *&
                        exp(xi) / ((exp(xi) + 1.0_prec)**2.0_prec)*&
                        (12.0_prec*log(2.0_prec*y) + 7.0_prec) 
        end function F_n_high_T 


        function F_rho_high_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi

                xi = phi - y

                res = Fermi(3, xi) + 3.0_prec * y * Fermi(2, xi) +&
                        2.5_prec*y*y * Fermi(1, xi) + .5_prec * y**3.0_prec * Fermi(0, xi)
                res = res / (y**3)

                res = res + exp(xi) / (exp(xi) + 1.0_prec) / 8.0_prec *&
                        (log(2.0_prec*y) + .25_prec) 
        end function F_rho_high_T 


        function F_p_high_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi

                xi = phi - y

                res = Fermi(3, xi) + 3.0_prec * y * Fermi(2, xi) +&
                        1.5_prec*y*y * Fermi(1, xi) - .5_prec * y**3.0_prec * Fermi(0, xi)
                res = res / (y**3)

                res = res + exp(xi) / (exp(xi) + 1.0_prec) / 8.0_prec *&
                        (3.0_prec * log(2.0_prec*y) + 1.75_prec) 
        end function F_p_high_T 


!-----------------------------
! Distribution functions for low temperature regime
! y -> +inf
!-----------------------------

        function F_n_low_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, zf, gf

                xi = phi - y
                gf = phi/y
                zf = sqrt(gf**2.0_prec - 1.0_prec)

                res = zf**3.0_prec / 3.0_prec + &
                        gint (f, 0.0_prec, xi)

                print*, gint(f, 0.0_prec, xi)


                contains

                function f (u)
                        real(prec) :: u, f
                        real(prec) :: gp, zp, gm, zm

                        gp = gf + u/y
                        zp = sqrt(gp**2.0_prec - 1.0_prec)

                        gm = gf - u/y
                        zm = sqrt(gm**2.0_prec - 1.0_prec)

                        f =(gp*zp - gm*zm) / (exp(u) + 1.0_prec)/y
                end function f

        end function F_n_low_T 


        function F_rho_low_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, zf, gf

                xi = phi - y
                gf = phi/y
                zf = sqrt(gf**2.0_prec - 1.0_prec)

                res = gf*2.0_prec * &
                        gint (f, 0.0_prec, xi)/y**2.0_prec

                res = res + (zf*gf * (2.0_prec *zf*zf + 1.0_prec) - &
                        log(zf + sqrt(zf*zf + 1.0_prec)))/8.0_prec



                contains

                function f (u)
                        real(prec) :: u, f
                        real(prec) :: gp, zp, gm, zm

                        gp = gf + u/y
                        zp = sqrt(gp**2.0_prec - 1.0_prec)

                        gm = gf - u/y
                        zm = sqrt(gm**2.0_prec - 1.0_prec)

                        f = u / (exp(u) + 1.0_prec) * &
                                (2.0_prec * (gf*gf + u*u/y/y)/ (zp + zm) - zp - zm)
                end function f
        end function F_rho_low_T 



        function F_p_low_T (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, zf, gf

                xi = phi - y
                gf = phi/y
                zf = sqrt(gf**2.0_prec - 1.0_prec)

                res = gf*4.0_prec * &
                        gint (f, 0.0_prec, xi)/y**2.0_prec

                res = res + (zf*gf * (2.0_prec *zf*zf - 3.0_prec) + &
                        3.0_prec*log(zf + sqrt(zf*zf + 1.0_prec)))/8.0_prec
                        

                contains

                function f (u)
                        real(prec) :: u, f
                        real(prec) :: gp, zp, gm, zm

                        gp = gf + u/y
                        zp = sqrt(gp**2.0_prec - 1.0_prec)

                        gm = gf - u/y
                        zm = sqrt(gm**2.0_prec - 1.0_prec)

                        f = u / (exp(u) + 1.0_prec) * &
                                (zp*zp + zm*zm + zp*zm)/(zp + zm)
                end function f

        end function F_p_low_T 


!-----------------------------
! Distribution functions without simplification
!-----------------------------

        function F_n (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, step


                xi = phi - y

                if (xi < 1.0_prec) then
                        step = 1.0_prec
                else
                        step = xi/10.0_prec
                endif
                
                if (y < 1.0_prec) then
                        res = gint_to_inf(f, 0.0_prec, step, 5)
                else
                        res = gint_to_inf(g, 0.0_prec, step, 5)/ (y**3.0_prec)
                endif


                contains

                function f (x)
                        real(prec) :: x, f

                        f = sqrt(x* (2.0_prec*y + x)) * (x+y) /&
                                (exp(x - xi) + 1.0_prec) / (y**3.0_prec)
                end function f

                function g (x)
                        real(prec) :: x, g

                        g = sqrt(x* (2.0_prec*y + x)) * (x+y) /&
                                (exp(x - xi) + 1.0_prec) 
                end function g
        end function F_n

        function F_rho (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, step


                xi = phi - y

                if (xi < 1.0_prec) then
                        step = 1.0_prec
                else
                        step = xi/10.0_prec
                endif
                
                if (y < 1.0_prec) then
                        res = gint_to_inf(f, 0.0_prec, step, 5)
                else
                        res = gint_to_inf(g, 0.0_prec, step, 5)/ (y**3.0_prec)
                endif

                contains

                function f (x)
                        real(prec) :: x, f

                        f = sqrt(x* (2.0_prec*y + x)) * ((x+y)**2.0_prec) /&
                                (exp(x - xi) + 1.0_prec) / (y**3.0_prec)
                end function f

                function g (x)
                        real(prec) :: x, g

                        g = sqrt(x* (2.0_prec*y + x)) * ((x+y)**2.0_prec) /&
                                (exp(x - xi) + 1.0_prec) 
                end function g
        end function F_rho


        function F_p (y, phi) result (res)
                implicit none

                real(prec) :: y, phi
                real(prec) :: res

                real(prec) :: xi, step


                xi = phi - y

                if (xi < 1.0_prec) then
                        step = 1.0_prec
                else
                        step = xi/10.0_prec
                endif

                if (y < 1.0_prec) then
                        res = gint_to_inf(f, 0.0_prec, step, 5)
                else
                        res = gint_to_inf(g, 0.0_prec, step, 5)/ (y**3.0_prec)
                endif

                contains

                function f (x)
                        real(prec) :: x, f

                        f = (x* (2.0_prec*y + x))**1.5_prec /&
                                (exp(x - xi) + 1.0_prec)/(y**3.0_prec)
                end function f

                function g (x)
                        real(prec) :: x, g

                        g = (x* (2.0_prec*y + x))**1.5_prec /&
                                (exp(x - xi) + 1.0_prec)
                end function g
        end function F_p


end module distributions
