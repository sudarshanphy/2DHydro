subroutine read_par()
        
        ! this subroutine reads from the parameter file and 
        ! passes the values as output
        use sim_data 
        implicit none

        ! local varaibles to read the parameter file.
        character(len=100) :: buffer, label
        integer :: pos
        integer, parameter :: no = 15 ! file unit number
        integer :: ios = 0
        integer :: line = 0 ! tracks number of lines in the paramter file

        ! set some default values, incase the user doesn't provide
        xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0
        nx = 10; ny = 10; restart_no = 0; restart_step = 0
        t0 = 0.0; tf = 0.5; dt = 1.0e-4; out_dt = 1.0e-3
        gamma = 1.4; cfl = 0.4; grav = 0.0
        xbctype = "periodic"; ybctype = "periodic"; problem = "none" 
        basenm = "def"; restart = .false.; usegrav = .false.;
        recon_method = "weno3"; flux_solver = "hllc"
 
        open(no, file='./par_input.par')

        ! ios is negative if an end of record condition is encountered or if
        ! an endfile condition was detected.  It is positive if an error was
        ! detected.  ios is zero otherwise.

        do while (ios == 0)
           read(no, '(A)', iostat=ios) buffer ! read the line
           if (ios == 0) then
              line = line + 1 ! count the line

              ! Find the first instance of whitespace.  Split label and data.
              pos = scan(buffer, '=')
              label = buffer(1:pos-2) ! this makes sure that we just read the label
              buffer = buffer(pos+1:) ! this is to read the value of the label

              select case (label) ! based on the label, now assign values to predefined varaibles
              case ('xmin')
                 read(buffer, *, iostat=ios) xmin
              case ('xmax')
                 read(buffer, *, iostat=ios) xmax
              case ('ymin')
                 read(buffer, *, iostat=ios) ymin
              case ('ymax')
                 read(buffer, *, iostat=ios) ymax
              case ('nx')
                 read(buffer, *, iostat=ios) nx
              case ('ny')
                 read(buffer, *, iostat=ios) ny
              case ('t0')
                 read(buffer, *, iostat=ios) t0
              case ('tf')
                 read(buffer, *, iostat=ios) tf
              case ('dt')
                 read(buffer, *, iostat=ios) dt
              case ('out_dt')
                 read(buffer, *, iostat=ios) out_dt
              case ('gamma')
                 read(buffer, *, iostat=ios) gamma
              case ('cfl')
                 read(buffer, *, iostat=ios) cfl
              case ('grav')
                 read(buffer, *, iostat=ios) grav
              case ('usegrav')
                 read(buffer, *, iostat=ios) usegrav
              case ('xbctype')
                 read(buffer, *, iostat=ios) xbctype
              case ('ybctype')
                 read(buffer, *, iostat=ios) ybctype
              case ('problem')
                 read(buffer, *, iostat=ios) problem
              case ('basenm')
                 read(buffer, *, iostat=ios) basenm
              case ('recon_method')
                 read(buffer, *, iostat=ios) recon_method
              case ('flux_solver')
                 read(buffer, *, iostat=ios) flux_solver
              case ('restart')
                 read(buffer, *, iostat=ios) restart
              case ('restart_no')
                 read(buffer, *, iostat=ios) restart_no
              case ('restart_step')
                 read(buffer, *, iostat=ios) restart_step
              case default
                ! for comments on the paramter file                
              end select
           end if
        end do
        close(no) !close the file
end subroutine
