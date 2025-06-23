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
       print *, "Sedov problem selected"
    end if
    do j = jlo, jhi
      do i = ilo, ihi
         ! initialize the fields
         solnVar(DENS_VAR, i,j) = 1.0e0
         solnVar(VELX_VAR, i,j) = 0.0e0
         solnVar(VELY_VAR, i,j) = 0.0e0
         solnVar(PRES_VAR, i,j) = 1.0e-1

         ! add a region with very high pressure
         if (sqrt(x(i)**2 + y(j)**2) <= 0.1) then
           solnVar(PRES_VAR, i,j) = 1.0e1
         end if
         call eos_gete(solnVar(:,i,j))
      end do
    end do
    deallocate(x,y)  
    nullify(solnVar)  
  end subroutine init_problem
end module sim_init
