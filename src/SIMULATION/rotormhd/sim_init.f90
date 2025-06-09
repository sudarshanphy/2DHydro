module sim_init
#include "param.h"
  implicit none
contains
  subroutine init_problem()

    use sim_data, only: problem, ihi, ilo, jhi, &
                        jlo, PI, grav, usegrav, gamma, &
                        mainVar, myrank
    use grid_func, only: get_coords
    use misc_module, only: to_upper
    use eos_module, only: eos_gete, eos_getp

    implicit none
    real, pointer :: solnVar(:,:,:)
    integer :: i, j, iInt, jInt
    real :: r
    real, allocatable :: x(:), y(:)

    allocate(x(ilo:ihi))
    allocate(y(jlo:jhi))

    x = 0.0
    y = 0.0
    
    call get_coords('x',ilo,ihi,x)
    call get_coords('y',jlo,jhi,y)

    solnVar(1:,ilo:,jlo:) => mainVar(1:,ilo:ihi,jlo:jhi)    
    solnVar(:,:,:) = 0.0
    
    ! setup take from "Limiting and divergence cleaning for continuous finite
    ! element discretizations of the MHD equations" by Dmitri Kuzmin
    ! page no 18
    if (myrank == MASTER_PROC) then 
       print *, "2D ROTOR MHD Problem"
    end if
    do j = jlo, jhi
      do i = ilo, ihi
         r = sqrt((x(i) - 0.5)**2 + (y(j) - 0.5)**2)
         if (r < 0.1) then
           solnVar(DENS_VAR, i,j) = 10.0
           solnVar(VELX_VAR, i,j) = 10.0 * (0.5 - y(j)) 
           solnVar(VELY_VAR, i,j) = 10.0 * (x(i) - 0.5)
         else if (r > 0.115) then
           solnVar(DENS_VAR, i,j) = 1.0
           solnVar(VELX_VAR, i,j) = 0.0
           solnVar(VELY_VAR, i,j) = 0.0
         else
           solnVar(DENS_VAR, i,j) = 1.0 + 9.0 * (0.115 - r)/(r - 0.1)
           solnVar(VELX_VAR, i,j) =  100.0 * (0.115 - r)/(r - 0.1) * (0.5 - y(j)) &
                                     / solnVar(DENS_VAR, i,j)
           solnVar(VELY_VAR, i,j) =  100.0 * (0.115 - r)/(r - 0.1) * (x(i) - 0.5) &
                                     / solnVar(DENS_VAR, i,j)
           ! need this to make sure density doen't blow up at initialization
           solnVar(DENS_VAR, i,j) = min(solnVar(DENS_VAR, i,j), 10.0)
         end if
         solnVar(PRES_VAR, i,j) = 0.5
#ifdef MHD                   
         solnVar(BMFX_VAR, i,j) = 2.5 / sqrt(4.0 * PI)
         solnVar(BMFY_VAR, i,j) = 0.0
#endif
         call eos_gete(solnVar(:,i,j)) 
      end do
    end do
    deallocate(x,y)  
    nullify(solnVar)  
  end subroutine init_problem
end module sim_init
