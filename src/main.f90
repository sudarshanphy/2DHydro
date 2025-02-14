program hydro
#include "header.h"
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
  real, allocatable, dimension(:, :, :) :: solnVar

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
  allocate(solnVar(xTpts, yTpts, NVAR_NUMBER))

  ! check if it is a restart or not
  if (.not. restart) then
    step = 0
    outputno = 0

    !initialize the problem
    call init_problem(xval, yval, solnVar(:,:,:))
                                             
    call write_output(t0, step, xval, yval, solnVar(ilo:ihi, jlo:jhi, :), &
                      outputno) 

  else

    ! restart a problem from an output file
    call restart_problem(restart_no, t0, solnVar)
    step = restart_step
    outputno = restart_no
    call write_output(t0, step, xval, yval, solnVar(ilo:ihi, jlo:jhi, :), &
                      outputno, .true.) 

  end if

  call applyBC_all(solnVar)
   
  time = t0
  timeio = time + out_dt
   
  ! evolution loop
  do while (time < tf)
    io_output = .false.
    
    ! compute dt
    dt = get_dt(solnVar(:,:,:))

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
    call glm(solnVar(:,:,BPSI_VAR), dt)
#endif 
    ! do a SSP RK2 step
    call RK2_SSP(dt, solnVar)
    ! apply BC
    call applyBC_all(solnVar)

    time = time + dt
    step = step + 1

    if (io_output) then
      outputno = outputno + 1
      call write_output(time, step, xval, yval, solnVar, outputno)
      print *, "$$$ Step: ", step," || time = ", time, " || output # = ", outputno," $$$" 

    end if
#ifdef MHD
    call glm(solnVar(:,:,BPSI_VAR), dt)
#endif
    print *, "Step: ", step," --> time = ", time , ", dt = ", dt
  end do
  
  deallocate(solnVar)
  deallocate(xval, yval)
end program hydro
