module config
        use :: precision
        implicit none 

        character(len = 100) :: output_dir
        character(len = 100) :: regime
        real(prec) :: Tmin, Tmax, mu, m, s, nu
        integer :: discr

        namelist /config_nml/ output_dir, Tmin, Tmax, discr, regime, mu, m, s, nu


contains

subroutine read_config (config_file)
        character (len = *) :: config_file

        ! fu - file unit
        ! io, if_exists - error flags
        integer :: fu, io
        logical :: if_exists


        ! check whether file exists
        inquire (file = config_file, exist = if_exists)

        if (.not. if_exists) then
                write (*, '("Error: configuration file ", a, " does not exist")') config_file
                close (fu)

                STOP 1
        end if

        ! open and read Namelist file
        open (action = 'read', file = config_file, iostat = io, newunit = fu)
        read (nml = config_nml, iostat = io, unit = fu)

        if (io /= 0) then
                write (*, '("Error: invalid namelist format")')
                close (fu)

                STOP 1
        end if
        
        close (fu)
end subroutine read_config


end module config
