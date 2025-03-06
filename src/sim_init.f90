module sim_init
#include "param.h"
  implicit none
contains
  subroutine init_problem()

    use sim_data, only: problem, Gpts, ihi, ilo, jhi, &
                        jlo, PI, grav, usegrav, gamma, &
                        mainVar
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
    y = 1.0
    
    call get_coords('x',ilo,ihi,x)
    call get_coords('y',jlo,jhi,y)

    solnVar(ilo:,jlo:,1:) => mainVar(ilo:ihi,jlo:jhi,1:)    
    solnVar(:,:,:) = 0.0

    select case (to_upper(trim(problem)))

    case ("SEDOV")
      print *, "Sedov problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 

           ! initialize the fields
           solnVar(i,j, DENS_VAR) = 1.0e0
           solnVar(i,j, VELX_VAR) = 0.0e0
           solnVar(i,j, VELY_VAR) = 0.0e0
           solnVar(i,j, PRES_VAR) = 1.0e-1

           ! add a region with very high pressure
           if (sqrt(x(iInt)**2 + y(jInt)**2) <= 0.1) then
             solnVar(i,j,PRES_VAR) = 1.0e1
           end if
           call eos_gete(solnVar(i,j,:))
        end do
      end do
    
    case ("CONSTANT")
      print *, "Constant problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 

           ! initialize the fields
           solnVar(i,j, DENS_VAR) = 1.0e0
           solnVar(i,j, VELX_VAR) = 0.0e0
           solnVar(i,j, VELY_VAR) = 0.0e0
           solnVar(i,j, PRES_VAR) = 1.0e0

           !solnVar(i,j,ENER_VAR) = eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END))
           call eos_gete(solnVar(i,j,:))
        end do
      end do

    case ("KH")
      print *, "Kelvin-Helmholtz (KH) problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 

           ! initialize the fields
           if (abs(y(jInt) - 0.5) > 0.25) then
             solnVar(i,j, DENS_VAR) = 1.0e0
             solnVar(i,j, VELX_VAR) = -0.5e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 2.5e0
           else
             solnVar(i,j, DENS_VAR) = 2.0e0
             solnVar(i,j, VELX_VAR) = 0.5e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 2.5e0
           end if
           ! give velocity perturbations
           if (abs(y(jInt) - 0.25) < 0.1) then
              solnVar(i,j, VELY_VAR) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
           if (abs(y(jInt) - 0.75) < 0.1) then
              solnVar(i,j, VELY_VAR) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
           call eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END)) 
        end do
      end do

    case ("RT")
      print *, "Rayleigh-Taylor (RT) problem selected"
      if (.not. usegrav) then
        print *, "Gravity is needed. Setting usegrav = T"
        usegrav = .true.
      end if
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 

           ! initialize the fields
           if (y(jInt) > 0.0) then
             solnVar(i,j, DENS_VAR) = 2.0e0
             solnVar(i,j, VELX_VAR) = -0.0e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 2.5e0 + grav * solnVar(i,j, DENS_VAR) * y(jInt)
           else
             solnVar(i,j, DENS_VAR) = 1.0e0
             solnVar(i,j, VELX_VAR) = 0.0e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 2.5e0 + grav * solnVar(i,j, DENS_VAR) * y(jInt)
           end if

           ! give velocity perturbation at the interface
           solnVar(i,j,VELY_VAR) = 0.01 * (1.0 + cos(4.0 * PI * x(iInt))) &
                               * (1.0 + cos(3.0 * PI * y(jInt))) / 4.0
                             
           call eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END)) 
        end do
      end do

    case ("RIEMANN")
      print *, "2D Riemann problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 

           ! initialize the fields
           if ((x(iInt) <= 0.8) .and. (y(jInt) <= 0.8)) then
             solnVar(i,j, DENS_VAR) = 0.138e0
             solnVar(i,j, VELX_VAR) = 1.206e0
             solnVar(i,j, VELY_VAR) = 1.206e0
             solnVar(i,j, PRES_VAR) = 0.029e0
           else if ((x(iInt) <= 0.8) .and. (y(jInt) > 0.8)) then
             solnVar(i,j, DENS_VAR) = 0.5323e0
             solnVar(i,j, VELX_VAR) = 1.206e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 0.3e0
           else if ((x(iInt) > 0.8) .and. (y(jInt) <= 0.8)) then
             solnVar(i,j, DENS_VAR) = 0.5323e0
             solnVar(i,j, VELX_VAR) = 0.0e0
             solnVar(i,j, VELY_VAR) = 1.206e0
             solnVar(i,j, PRES_VAR) = 0.3e0
           else
             solnVar(i,j, DENS_VAR) = 1.5e0
             solnVar(i,j, VELX_VAR) = 0.0e0
             solnVar(i,j, VELY_VAR) = 0.0e0
             solnVar(i,j, PRES_VAR) = 1.5e0
           end if
           call eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END)) 
        end do
      end do

    case ("OTMHD")
      print *, "2D Orszag-Tang MHD Vortex"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 
           solnVar(i,j, DENS_VAR) = gamma**2
           solnVar(i,j, VELX_VAR) = -sin(y(j))
           solnVar(i,j, VELY_VAR) = sin(x(i))
           solnVar(i,j, PRES_VAR) = gamma
#ifdef MHD
           solnVar(i,j, BMFX_VAR) = -sin(y(j))
           solnVar(i,j, BMFY_VAR) = sin(2.0*x(i))
#endif 
           call eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END)) 
        end do
      end do

    case ("ROTORMHD")
      print *, "2D ROTOR MHD Problem"
      ! setup take from "Limiting and divergence cleaning for continuous finite
      ! element discretizations of the MHD equations" by Dmitri Kuzmin
      ! page no 18
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i 
           jInt = j 
           r = sqrt((x(iInt) - 0.5)**2 + (y(jInt) - 0.5)**2)
           if (r < 0.1) then
             solnVar(i,j, DENS_VAR) = 10.0
             solnVar(i,j, VELX_VAR) = 10.0 * (0.5 - y(jInt)) 
             solnVar(i,j, VELY_VAR) = 10.0 * (x(iInt) - 0.5)
           else if (r > 0.115) then
             solnVar(i,j, DENS_VAR) = 1.0
             solnVar(i,j, VELX_VAR) = 0.0
             solnVar(i,j, VELY_VAR) = 0.0
           else
             solnVar(i,j, DENS_VAR) = 1.0 + 9.0 * (0.115 - r)/(r - 0.1)
             solnVar(i,j, VELX_VAR) =  100.0 * (0.115 - r)/(r - 0.1) * (0.5 - y(jInt)) &
                                       / solnVar(i,j, DENS_VAR)
             solnVar(i,j, VELY_VAR) =  100.0 * (0.115 - r)/(r - 0.1) * (x(iInt) - 0.5) &
                                       / solnVar(i,j, DENS_VAR)
             ! need this to make sure density doen't blow up at initialization
             solnVar(i,j,DENS_VAR) = min(solnVar(i,j,DENS_VAR), 10.0)
           end if
           solnVar(i,j, PRES_VAR) = 0.5
#ifdef MHD
           solnVar(i,j, BMFX_VAR) = 2.5 / sqrt(4.0 * PI)
           solnVar(i,j, BMFY_VAR) = 0.0
#endif
           call eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END)) 
        end do
      end do
    case default
      print *, "such problem not defined"
    end select

  deallocate(x,y)  
  nullify(solnVar)  
  end subroutine init_problem
end module sim_init
