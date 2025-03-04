module sim_data
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

  ! array which has the solution
  real, save, dimension(:,:,:), target, allocatable :: mainVar

  ! only used by MHD runs 
  real, save :: ch !speed for divergence cleaning
end module sim_data
