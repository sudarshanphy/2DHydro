module sim_data
  use mpi
  implicit none
  real, save :: xmin, xmax, ymin, ymax, dx, dy
  real, save :: t0, tf, dt, out_dt
  real, save :: gamma, cfl, grav
  integer, save :: nx, ny, Gpts, xTpts, yTpts, restart_no, &
                   restart_step 
  character(len=256), save :: xbctype, ybctype, problem, &
                              basenm, recon_method, flux_solver
  logical, save :: restart, usegrav
  real, parameter :: PI = 4.0 * atan(1.0)
  real, parameter :: smallf = 1.0e-30

  ! array which has the solution
  real, save, dimension(:,:,:), target, allocatable :: mainVar

  ! only used by MHD runs 
  real, save :: ch !speed for divergence cleaning


  ! MPI stuff
  integer, save :: myrank, nprocs, comm, ierr, errcode
  integer, save :: xrank, yrank
  integer, dimension(MPI_STATUS_SIZE) :: status1 

  integer, save :: ilo, ihi, jlo, jhi, iGlo, iGhi, jGlo, jGhi
  integer, save :: xblk, yblk, lnx, lny, lxTpts, lyTpts
  ! local block xmin and xmax
  real, save :: lxmin, lxmax, lymin, lymax
  ! blocks to the left, right, top, bottom
  integer, save :: l_blk, r_blk, t_blk, b_blk
  logical, save :: at_xboundary, at_yboundary
end module sim_data
