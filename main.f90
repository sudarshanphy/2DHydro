program hydro
#include "header.h"

  use sim_data
  use sim_init, only: init_problem
  use sim_restart, only: restart_problem
  use io_module, only: write_output
  use applyBC_module, only: applyBC_all
  use rk2_module, only: RK2_SSP
  use misc_module, only: get_dt

  implicit none

  real, allocatable, dimension(:) :: xval, yval
  real, allocatable, dimension(:, :) :: dens, velx, vely, pres, ener

#ifdef MHD
  real, allocatable, dimension(:, :) :: bmfx, bmfy, bpsi
#endif   

  real :: time, dtime, timeio
  integer ::  step, outputno
  logical :: io_output
  
  !read the parameter file
  call read_par()
  
  !allocate space for the interior grid points
  allocate(xval(nx), yval(ny))

  !initialize the grid
  call grid_init(xval, yval)

  !allocate the fields
  allocate(dens(xTpts, yTpts), velx(xTpts, yTpts), &
           vely(xTpts, yTpts), pres(xTpts, yTpts), ener(xTpts, yTpts))
#ifdef MHD  
  allocate(bmfx(xTpts, yTpts), bmfy(xTpts, yTpts), &
           bpsi(xTpts, yTpts))
#endif

  ! check if it is a restart or not
  if (.not. restart) then
    step = 0
    outputno = 0

    !initialize the problem
#ifdef MHD
    call init_problem(xval, yval, dens, velx, vely, pres, ener, &
                                              bmfx, bmfy, bpsi )
#else
    call init_problem(xval, yval, dens, velx, vely, pres, ener)
#endif
                                             
    call write_output(t0, step, xval, yval, &
                      dens(ilo:ihi, jlo:jhi), velx(ilo:ihi, jlo:jhi), &
                      vely(ilo:ihi, jlo:jhi), pres(ilo:ihi, jlo:jhi), &
                      ener(ilo:ihi, jlo:jhi), &
#ifdef MHD
                      bmfx(ilo:ihi, jlo:jhi), bmfy(ilo:ihi, jlo:jhi), &
                      bpsi(ilo:ihi, jlo:jhi), &
#endif
                      outputno) 

  else

    ! restart a problem from an output file
#ifdef MHD
    call restart_problem(restart_no, t0, dens, velx, vely, pres, ener, &
                                                     bmfx, bmfy, bpsi )
#else
    call restart_problem(restart_no, t0, dens, velx, vely, pres, ener)
#endif

    step = restart_step
    outputno = restart_no
    call write_output(t0, step, xval, yval, &
                      dens(ilo:ihi, jlo:jhi), velx(ilo:ihi, jlo:jhi), &
                      vely(ilo:ihi, jlo:jhi), pres(ilo:ihi, jlo:jhi), &
                      ener(ilo:ihi, jlo:jhi),                         &
#ifdef MHD
                      bmfx(ilo:ihi, jlo:jhi), bmfy(ilo:ihi, jlo:jhi), &
                      bpsi(ilo:ihi, jlo:jhi), &
#endif
                      outputno, .true.) 

  end if

#ifdef MHD
  call applyBC_all(dens, velx, vely, pres, ener, bmfx, bmfy, bpsi)
#else
  call applyBC_all(dens, velx, vely, pres, ener)
#endif  
   
  time = t0
  timeio = time + out_dt
   
  ! evolution loop
  do while (time < tf)
    io_output = .false.
    
    ! compute dt
#ifdef MHD
    dt = get_dt(dens, velx, vely, pres, bmfx, bmfy)
#else    
    dt = get_dt(dens, velx, vely, pres) 
#endif

    if (time + dt > tf) then
      dt = tf - time
      io_output = .true.
    end if
    
    if (time + dt > timeio) then
      dt = timeio - time
      timeio = timeio + out_dt
      io_output = .true.
    end if


#ifdef MHD
    ! glm update psi 
    call glm(bpsi, dt)
    ! do a SSP RK2 step
    call RK2_SSP(dt, dens, velx, vely, pres, ener, bmfx, bmfy, bpsi)
    ! apply BC
    call applyBC_all(dens, velx, vely, pres, ener, bmfx, bmfy, bpsi)
#else
    call RK2_SSP(dt, dens, velx, vely, pres, ener)
    call applyBC_all(dens, velx, vely, pres, ener) 
#endif  
    time = time + dt
    step = step + 1

    if (io_output) then
      outputno = outputno + 1
      call write_output(time, step, xval, yval, &
                        dens(ilo:ihi, jlo:jhi), velx(ilo:ihi, jlo:jhi), &
                        vely(ilo:ihi, jlo:jhi), pres(ilo:ihi, jlo:jhi), &
                        ener(ilo:ihi, jlo:jhi),                         & 
#ifdef MHD
                        bmfx(ilo:ihi, jlo:jhi), bmfy(ilo:ihi, jlo:jhi), &
                        bpsi(ilo:ihi, jlo:jhi), &
#endif
                                                                  outputno)
      print *, "$$$ Step: ", step," || time = ", time, " || output # = ", outputno," $$$" 

    end if
#ifdef MHD
    call glm(bpsi, dt)
#endif
    print *, "Step: ", step," --> time = ", time , ", dt = ", dt
  end do
  
  deallocate(dens, velx, vely, pres, ener)
#ifdef MHD
  deallocate(bmfx, bmfy, bpsi)
#endif
  deallocate(xval, yval)
end program hydro
