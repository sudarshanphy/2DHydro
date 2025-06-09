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
       print *, "Rayleigh-Taylor (RT) problem selected"
    end if
    if (.not. usegrav) then
      usegrav = .true.
      if (myrank == MASTER_PROC) then
         print *, "Gravity is needed. Setting usegrav = T"
      end if
    end if
    do j = jlo, jhi
      do i = ilo, ihi
         ! initialize the fields
         if (y(j) > 0.0) then
           solnVar(DENS_VAR, i,j) = 2.0e0
           solnVar(VELX_VAR, i,j) = -0.0e0
           solnVar(VELY_VAR, i,j) = 0.0e0
           solnVar(PRES_VAR, i,j) = 2.5e0 + grav * solnVar(DENS_VAR, i,j) * y(j)
         else
           solnVar(DENS_VAR, i,j) = 1.0e0
           solnVar(VELX_VAR, i,j) = 0.0e0
           solnVar(VELY_VAR, i,j) = 0.0e0
           solnVar(PRES_VAR, i,j) = 2.5e0 + grav * solnVar(DENS_VAR, i,j) * y(j)
         end if

         ! give velocity perturbation at the interface
         solnVar(VELY_VAR, i,j) = 0.01 * (1.0 + cos(4.0 * PI * x(i))) &
                             * (1.0 + cos(3.0 * PI * y(j))) / 4.0
                           
         call eos_gete(solnVar(:,i,j)) 
      end do
    end do
    deallocate(x,y)  
    nullify(solnVar)  
  end subroutine init_problem
end module sim_init
