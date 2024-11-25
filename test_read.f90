program test_read
  use sim_data
  implicit none
  
  call read_par()
  print *, "problem = ", trim(problem)
  print *, "basenm = ", trim(basenm)
  print *, "xbctype, ybctype = ", trim(xbctype), ", ", trim(ybctype)
  print *, "xmin, xmax, ymin, ymax = ", xmin, ", ", xmax, ", ", &
                                        ymin, ", ", ymax
  print *, "nx, ny = ", nx, ", ", ny
  print *, "t0, tf, dt, out_dt = ", t0, ", ", tf, ", ", dt, ", ", out_dt                                    
  print *, "gamma, cfl = ", gamma, ", ", cfl
  print *, "restart, retsart_no, restart_step = ", restart, ", ", &
                     restart_no, ", ", restart_step

end program test_read
