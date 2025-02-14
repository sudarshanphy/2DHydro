module misc_module
#include "param.h"
  implicit none
contains

    function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
         implicit none
    
         character(len=*), intent(in) :: strIn
         character(len=len(strIn)) :: strOut
         integer :: i,j
    
         do i = 1, len(strIn)
              j = iachar(strIn(i:i))
              if (j>= iachar("a") .and. j<=iachar("z") ) then
                   strOut(i:i) = achar(iachar(strIn(i:i))-32)
              else
                   strOut(i:i) = strIn(i:i)
              end if
         end do
    
    end function to_upper

    function get_dt(solnVar) result(dt)
      use sim_data, only: gamma, xTpts, yTpts, dx, dy, cfl
#ifdef MHD
      use sim_data, only: ch
#endif
      implicit none
      real, dimension(xTpts, yTpts, NVAR_NUMBER) :: solnVar
      real, dimension(xTpts, yTpts) :: xsmax, ysmax
      real :: cs, xcmax, ycmax   !sound speed
#ifdef MHD  
      real :: cax, cay, cfx, cfy, B2, cB2 !alfven wave and magnetosonic wave
#endif
      real :: dt
      integer :: i, j
      real :: max_xsmax, max_ysmax 
      
      do j = 1, yTpts
        do i = 1, xTpts
           cs = sqrt(gamma * solnVar(i,j,PRES_VAR) / solnVar(i,j,DENS_VAR))

           if (cs < 0.0) then
             print *, "Imaginary sound speed!"
             print *, "At i,j, with gamma, pres, dens = ", i, j, gamma, &
                           solnVar(i,j,PRES_VAR), solnVar(i,j,DENS_VAR)
           endif
           xcmax = cs
           ycmax = cs
#ifdef MHD
           cax = solnVar(i,j,BMFX_VAR) / sqrt(solnVar(i,j,DENS_VAR))
           cay = solnVar(i,j,BMFY_VAR) / sqrt(solnVar(i,j,DENS_VAR))
           B2 =  sum(solnVar(i,j,BMFX_VAR:BMFZ_VAR)*solnVar(i,j,BMFX_VAR:BMFZ_VAR))
           cB2 = B2 / solnVar(i,j,DENS_VAR)
           cfx = sqrt(0.5 * ((cs + cB2) + sqrt((cs + cB2)**2 - 4.0 * cs * cax * cax)))
           cfy = sqrt(0.5 * ((cs + cB2) + sqrt((cs + cB2)**2 - 4.0 * cs * cay * cay)))
           xcmax = cfx
           xcmax = cfy
#endif
           xsmax(i,j) = xcmax + abs(solnVar(i,j,VELX_VAR))
           ysmax(i,j) = ycmax + abs(solnVar(i,j,VELY_VAR))
        end do
      end do
       
      max_xsmax = maxval(xsmax); max_ysmax = maxval(ysmax)
      
#ifdef MHD
      !speed for divergence cleaning
      ch = max(max_xsmax, max_ysmax)
#endif

      !print *, "max xsmax, ysmax = ", max_xsmax, max_ysmax 
      dt = min(cfl * dx/max_xsmax, cfl * dy/max_ysmax) 

    end function get_dt

end module misc_module
