program hydro
#include "param.h"

  use sim_data
  use sim_init, only: init_problem
  use sim_restart, only: restart_problem
  use io_module, only: write_output
  use applyBC_module, only: applyBC_all
  use rk2_module, only: RK2_SSP
  use misc_module, only: get_dt 

  implicit none

  real, allocatable, dimension(:) :: xval, yval

  real :: time, dtime, timeio
  integer ::  step, outputno
  logical :: io_output
  
  real :: timer_start, timer_stop, time_per_step, timer_step0, &
          timer_step1
 
  real :: maxdivB

  !read the parameter file
  call read_par()
  
  !allocate space for the interior grid points
  allocate(xval(nx), yval(ny))

  !initialize the grid
  call grid_init(xval, yval)

  ! check if it is a restart or not
  if (.not. restart) then
    step = 0
    outputno = 0

    !initialize the problem
    call init_problem(xval, yval)
                                             
    call write_output(t0, step, xval, yval, outputno) 

  else

    ! restart a problem from an output file
    call restart_problem(restart_no, t0)
    step = restart_step
    outputno = restart_no
    call write_output(t0, step, xval, yval, &
                      outputno, .true.) 

  end if
  
  call CPU_TIME(timer_start)
  timer_step0 = 0
  call applyBC_all()
   
  time = t0
  timeio = time + out_dt
  ! evolution loop
  do while (time < tf)

    io_output = .false.
    
    ! compute dt
    call get_dt(dt)

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
    call glm(dt)
#endif 
    ! do a SSP RK2 step
    call RK2_SSP(dt)
    ! apply BC
    call applyBC_all()

    time = time + dt
    step = step + 1

    if (io_output) then

      call CPU_TIME(timer_stop)
      timer_step1 = step
      time_per_step = (timer_stop - timer_start) &
                     /(timer_step1 - timer_step0)

      print *, "# step0 = ", timer_step0, " step1 = ", timer_step1
      print *, "# time per step = ", time_per_step

      timer_step0 = step

      outputno = outputno + 1
      call write_output(time, step, xval, yval, outputno)
      print *, "$$$ Step: ", step," || time = ", time, " || output # = ", outputno," $$$" 
      call CPU_TIME(timer_start)
    end if
#ifdef MHD
    call glm(dt)
#endif
    print *, "Step: ", step," --> time = ", time , ", dt = ", dt
  end do
  
  deallocate(mainVar)
  deallocate(xval, yval)
end program hydro
