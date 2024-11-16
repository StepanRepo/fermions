program fermions
        use :: precision
        use :: config
        use :: gauss_int
        use :: distributions
        use :: sle_solve
implicit none

real(prec), dimension(:), allocatable :: T
real(prec), dimension(:, :), allocatable :: u, ul, res
integer i

integer :: fu, io
logical :: is_opened

call configure()

write(*, "('Input values:')")

write(*, *)""
write(*, "('T_min         : ', es16.2, ' K')") Tmin
write(*, "('T_max         : ', es16.2, ' K')") Tmax
write(*, "('Discretization: ', i16)") discr
write(*, *)""
write(*, "('Chem. potent. : ', es16.8, ' erg/g')") mu
write(*, *)""
write(*, "('Mass          : ', es16.8, ' g')") m
if (m == 0.0_prec) &
        write(*, "('Frequency     : ', es16.8, ' Hz')") nu
write(*, "('Spin          : ', f16.1)") s
write(*, *)""
write(*, "('Regime        : ', a16)") trim(regime)
write(*, "('Precision     : ', i16)") prec
write(*, "('Output dir    : ', a16)") trim(output_dir)
write(*, *)""


allocate(res(0:discr-1, 0:2), T(0:discr-1))
allocate(u(0:discr-1, 0:1), ul(0:discr-1, 0:1))



! Define discretization grid
do i = 0, discr-1
        T(i) = real(i, prec)*(Tmax-Tmin)/real(discr, prec) + Tmin
enddo

u(:, 0) = T
u(:, 1) = mu

! convert it to unitless coordinates
ul = unit2ul(u)


! Do the calculations
if (trim(regime) == "normal") then
        do i = 0, discr-1
                res(i, 0) = F_n(ul(i, 0), ul(i, 1))
                res(i, 1) = F_rho(ul(i, 0), ul(i, 1))
                res(i, 2) = F_p(ul(i, 0), ul(i, 1))
        enddo
else if (trim(regime) == "lowT") then
        do i = 0, discr-1
                res(i, 0) = F_n_low_T(ul(i, 0), ul(i, 1))
                res(i, 1) = F_rho_low_T(ul(i, 0), ul(i, 1))
                res(i, 2) = F_p_low_T(ul(i, 0), ul(i, 1))
        enddo
else if (trim(regime) == "highT") then
        do i = 0, discr-1
                res(i, 0) = F_n_high_T(ul(i, 0), ul(i, 1))
                res(i, 1) = F_rho_high_T(ul(i, 0), ul(i, 1))
                res(i, 2) = F_p_high_T(ul(i, 0), ul(i, 1))
        enddo
else
        write (*, "('Unknown regime: ', a)") trim(regime)
endif

! convert distribution values to physical ones
!res = distr2phys(res)


! open output file
fu = 1
open (action = "write", file = trim(output_dir)//"/"//trim(regime)//".dat",&
        iostat = io, unit = fu)

! if the program can't open desirable file
! it creates new temporary file to save results
inquire(unit = fu, opened = is_opened)

if (.not. is_opened) then
        write(*, "('Can`t open file: ', a)") trim(output_dir)//"/"//trim(regime)
        write(*, "('Saving to ./temp.dat')")

        close (fu)

        fu = 2
        open (action = "write", file = "temp.dat", iostat = io, unit = fu)
endif


! save results
do i = 0, discr-1
        write(fu, *) u(i, :), res(i, :)
enddo

close (fu)


contains


subroutine configure()
        implicit none

        character(len = 100):: arg, config_file = ""
        integer :: i, idx

        idx = 1
        do
                if (idx > command_argument_count()) exit
                call get_command_argument(idx, arg)

                select case (arg)

                case ("-h", "--help")
                        write(*, "(a)") "Command-line options:"
                        write(*, "(a)") "  -h, --help             print this help message"
                        write(*, "(a)") "  -c, --cfg {file_name}  set configuration file"

                        STOP 0

                case ("-c", "--cfg")
                        idx = idx+1
                        call get_command_argument(idx, arg)

                        config_file = arg

                case default
                        write(*, "(2a)") "Unknown command line option: ", trim(arg)
                end select

                idx = idx+1
        end do

        if (trim(config_file) /= "") then
                call read_config(trim(config_file)) ! read custom config
        else
                call read_config("default.cfg") ! read default config
        endif

end subroutine configure



end program 
