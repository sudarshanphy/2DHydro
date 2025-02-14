module sim_init
#include "header.h"
#include "param.h"
  implicit none
contains
  subroutine init_problem(x, y, solnVar)

    use sim_data, only: problem, xTpts, yTpts, &
                        nx, ny, Gpts, ihi, ilo, jhi, &
                        jlo, PI, grav, usegrav, gamma

    use misc_module, only: to_upper
    use eos_module, only: eos_gete
    implicit none
    real, dimension(nx), intent(in) :: x
    real, dimension(ny), intent(in) :: y
    real, dimension(xTpts, yTpts, NVAR_NUMBER), intent(inout) :: solnVar
    integer :: i, j, iInt, jInt
    real :: r
    
    solnVar(:,:,:) = 0.0

    select case (to_upper(trim(problem)))

    case ("SEDOV")
      print *, "Sedov problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           solnVar(i,j, DENS_VAR) = 1.0e0
           solnVar(i,j, VELX_VAR) = 0.0e0
           solnVar(i,j, VELY_VAR) = 0.0e0
           solnVar(i,j, PRES_VAR) = 1.0e-1

           ! add a region with very high pressure
           if (sqrt(x(iInt)**2 + y(jInt)**2) <= 0.1) then
             solnVar(i,j,PRES_VAR) = 1.0e1
           end if
           solnVar(i,j,ENER_VAR) = eos_gete(solnVar(i,j,NVAR_BEGIN:NVAR_END-1)) 
        end do
      end do
#if 0
    case ("KH")
      print *, "Kelvin-Helmholtz (KH) problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           if (abs(y(jInt) - 0.5) > 0.25) then
             dens(i,j) = 1.0e0
             velx(i,j) = -0.5e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0
           else
             dens(i,j) = 2.0e0
             velx(i,j) = 0.5e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0
           end if
#ifdef MHD
           bmfx(i,j) = 0.0e0
           bmfy(i,j) = 0.0e0
           bpsi(i,j) = 0.0e0
#endif 
           ! give velocity perturbations
           if (abs(y(jInt) - 0.25) < 0.1) then
              vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
           if (abs(y(jInt) - 0.75) < 0.1) then
              vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
#ifdef MHD
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/)) 
#else
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
#endif
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
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           if (y(jInt) > 0.0) then
             dens(i,j) = 2.0e0
             velx(i,j) = -0.0e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0 + grav * dens(i,j) * y(jInt)
           else
             dens(i,j) = 1.0e0
             velx(i,j) = 0.0e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0 + grav * dens(i,j) * y(jInt)
           end if
#ifdef MHD
           bmfx(i,j) = 0.0e0
           bmfy(i,j) = 0.0e0
           bpsi(i,j) = 0.0e0
#endif 
           ! give velocity perturbation at the interface
           !if (abs(y(jInt)) < 0.1) then
              !vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
              vely(i,j) = 0.01 * (1.0 + cos(4.0 * PI * x(iInt))) &
                               * (1.0 + cos(3.0 * PI * y(jInt))) / 4.0
           !end if
#ifdef MHD
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/)) 
#else
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
#endif
        end do
      end do

    case ("RIEMANN")
      print *, "2D Riemann problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           if ((x(iInt) <= 0.8) .and. (y(jInt) <= 0.8)) then
             dens(i,j) = 0.138e0
             velx(i,j) = 1.206e0
             vely(i,j) = 1.206e0
             pres(i,j) = 0.029e0
           else if ((x(iInt) <= 0.8) .and. (y(jInt) > 0.8)) then
             dens(i,j) = 0.5323e0
             velx(i,j) = 1.206e0
             vely(i,j) = 0.0e0
             pres(i,j) = 0.3e0
           else if ((x(iInt) > 0.8) .and. (y(jInt) <= 0.8)) then
             dens(i,j) = 0.5323e0
             velx(i,j) = 0.0e0
             vely(i,j) = 1.206e0
             pres(i,j) = 0.3e0
           else
             dens(i,j) = 1.5e0
             velx(i,j) = 0.0e0
             vely(i,j) = 0.0e0
             pres(i,j) = 1.5e0
           end if
#ifdef MHD
           bmfx(i,j) = 0.0e0
           bmfy(i,j) = 0.0e0
           bpsi(i,j) = 0.0e0
#endif 
#ifdef MHD
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/)) 
#else
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
#endif
        end do
      end do

    case ("OTMHD")
      print *, "2D Orszag-Tang MHD Vortex"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts
           dens(i,j) = gamma**2
           pres(i,j) = gamma
           velx(i,j) = -sin(y(j))
           vely(i,j) = sin(x(i))
#ifdef MHD
           bmfx(i,j) = -sin(y(j))
           bmfy(i,j) = sin(2.0*x(i))
           bpsi(i,j) = 0.0e0
#endif 
#ifdef MHD
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/)) 
#else
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
#endif
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
           iInt = i - Gpts
           jInt = j - Gpts
           r = sqrt((x(iInt) - 0.5)**2 + (y(jInt) - 0.5)**2)
           if (r < 0.1) then
             dens(i,j) = 10.0
             velx(i,j) = 10.0 * (0.5 - y(jInt)) 
             vely(i,j) = 10.0 * (x(iInt) - 0.5)
           else if (r > 0.115) then
             dens(i,j) = 1.0
             velx(i,j) = 0.0
             vely(i,j) = 0.0
           else
             dens(i,j) = 1.0 + 9.0 * (0.115 - r)/(r - 0.1)
             velx(i,j) =  100.0 * (0.115 - r)/(r - 0.1) * (0.5 - y(jInt)) / dens(i,j)
             vely(i,j) =  100.0 * (0.115 - r)/(r - 0.1) * (x(iInt) - 0.5) / dens(i,j)
             ! need this to make sure density doen't blow up at initialization
             dens(i,j) = min(dens(i,j), 10.0)
           end if
           pres(i,j) = 0.5
           bmfx(i,j) = 2.5 / sqrt(4.0 * PI)
           bmfy(i,j) = 0.0
           bpsi(i,j) = 0.0e0
#ifdef MHD
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/)) 
#else
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
#endif
        end do
      end do
#endif

    case default
      print *, "such problem not defined"
    end select    
     
  end subroutine init_problem
end module sim_init
