module sim_data
#include "header.h"

  implicit none
  real, save :: xmin, xmax, ymin, ymax, dx, dy
  real, save :: t0, tf, dt, out_dt
  real, save :: gamma, cfl, grav 
  integer, save :: nx, ny, ilo, ihi, jlo, jhi, Gpts, xTpts, yTpts, &
                   restart_no, restart_step
  character(len=256), save :: xbctype, ybctype, problem, &
                              basenm, recon_method, flux_solver
  logical, save :: restart, usegrav
  real, parameter :: PI = 4.0 * atan(1.0)
  real, parameter :: smallf = 1.0e-30 
#ifdef MHD
  real, save :: ch !speed for divergence cleaning
#endif
end module sim_data
