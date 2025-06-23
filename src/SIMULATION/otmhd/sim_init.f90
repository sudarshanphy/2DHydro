module sim_init
#include "param.h"
  implicit none
contains
  subroutine init_problem()

    use sim_data, only: ihi, ilo, jhi, &
                        jlo, PI, grav, usegrav, gamma, &
                        mainVar, myrank
    use grid_func, only: get_coords
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
    
    call get_coords(IAXIS,ilo,ihi,x)
    call get_coords(JAXIS,jlo,jhi,y)

    solnVar(1:,ilo:,jlo:) => mainVar(1:,ilo:ihi,jlo:jhi)    
    solnVar(:,:,:) = 0.0
    
    if (myrank == MASTER_PROC) then 
       print *, "2D Orszag-Tang(OT) MHD Vortex problem"
    end if
    do j = jlo, jhi
      do i = ilo, ihi
         solnVar(DENS_VAR, i,j) = gamma**2
         solnVar(VELX_VAR, i,j) = -sin(y(j))
         solnVar(VELY_VAR, i,j) = sin(x(i))
         solnVar(PRES_VAR, i,j) = gamma
#ifdef MHD
         solnVar(BMFX_VAR, i,j) = -sin(y(j))
         solnVar(BMFY_VAR, i,j) = sin(2.0*x(i))
#endif 
         call eos_gete(solnVar(:,i,j)) 
      end do
    end do
    deallocate(x,y)  
    nullify(solnVar)  
  end subroutine init_problem
end module sim_init
