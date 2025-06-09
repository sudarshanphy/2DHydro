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
    
    if (myrank == MASTER_PROC) then 
       print *, "2D Riemann problem selected"
    end if
    do j = jlo, jhi
      do i = ilo, ihi
         ! initialize the fields
         if ((x(i) <= 0.8) .and. (y(j) <= 0.8)) then
           solnVar(DENS_VAR, i,j) = 0.138e0
           solnVar(VELX_VAR, i,j) = 1.206e0
           solnVar(VELY_VAR, i,j) = 1.206e0
           solnVar(PRES_VAR, i,j) = 0.029e0
         else if ((x(i) <= 0.8) .and. (y(j) > 0.8)) then
           solnVar(DENS_VAR, i,j) = 0.5323e0
           solnVar(VELX_VAR, i,j) = 1.206e0
           solnVar(VELY_VAR, i,j) = 0.0e0
           solnVar(PRES_VAR, i,j) = 0.3e0
         else if ((x(i) > 0.8) .and. (y(j) <= 0.8)) then
           solnVar(DENS_VAR, i,j) = 0.5323e0
           solnVar(VELX_VAR, i,j) = 0.0e0
           solnVar(VELY_VAR, i,j) = 1.206e0
           solnVar(PRES_VAR, i,j) = 0.3e0
         else
           solnVar(DENS_VAR, i,j) = 1.5e0
           solnVar(VELX_VAR, i,j) = 0.0e0
           solnVar(VELY_VAR, i,j) = 0.0e0
           solnVar(PRES_VAR, i,j) = 1.5e0
         end if
         call eos_gete(solnVar(:,i,j)) 
      end do
    end do
    deallocate(x,y)  
    nullify(solnVar)  
  end subroutine init_problem
end module sim_init
