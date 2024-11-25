program hydro
  use sim_data
  use sim_init, only: init_problem
  use io_module, only: write_output
  use applyBC_module, only: applyBC_all
  use rk2_module, only: RK2_SSP
  use misc_module, only: get_dt 
  implicit none
  real, allocatable, dimension(:) :: xval, yval
  real, allocatable, dimension(:, :) :: dens, velx, vely, pres, ener
   
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
  if (.not. restart) then
    step = 0
    outputno = 0
    !initialize the problem
    call init_problem(xval, yval, dens, velx, vely, pres, ener)
    call write_output(t0, step, xval, yval, &
                      dens(ilo:ihi, jlo:jhi), velx(ilo:ihi, jlo:jhi), &
                      vely(ilo:ihi, jlo:jhi), pres(ilo:ihi, jlo:jhi), &
                      ener(ilo:ihi, jlo:jhi), outputno) 

  end if
  call applyBC_all(dens, velx, vely, pres, ener)
   
  time = t0
  timeio = time + out_dt
   
  ! evolution loop
  do while (time < tf)
    io_output = .false.
    
    ! compute dt
    dt = get_dt(dens, velx, vely, pres)

    if (time + dt > tf) then
      dt = tf - time
      io_output = .true.
    end if
    
    if (time + dt > timeio) then
      dt = timeio - time
      timeio = timeio + out_dt
      io_output = .true.
    end if

    call RK2_SSP(dens, velx, vely, pres, ener, dt)
    call applyBC_all(dens, velx, vely, pres, ener)
    time = time + dt
    step = step + 1

    if (io_output) then
      outputno = outputno + 1
      call write_output(time, step, xval, yval, &
                        dens(ilo:ihi, jlo:jhi), velx(ilo:ihi, jlo:jhi), &
                        vely(ilo:ihi, jlo:jhi), pres(ilo:ihi, jlo:jhi), &
                        ener(ilo:ihi, jlo:jhi), outputno)
      print *, "$$$ Step: ", step," || time = ", time, " || output # = ", outputno," $$$" 

    end if
    print *, "Step: ", step," --> time = ", time , ", dt = ", dt
  end do
  
  deallocate(dens, velx, vely, pres, ener)
  deallocate(xval, yval)
end program hydro
