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

    subroutine get_dt(dt)
      use sim_data, only: gamma, xTpts, yTpts, dx, dy, cfl, &
                          mainVar
#ifdef MHD
      use sim_data, only: ch
#endif
      implicit none
      real, dimension(xTpts, yTpts) :: xsmax, ysmax
      real :: cs, xcmax, ycmax   !sound speed
#ifdef MHD  
      real :: cax, cay, cfx, cfy, B2, cB2 !alfven wave and magnetosonic wave
#endif
      real, intent(out) :: dt
      integer :: i, j
      real :: max_xsmax, max_ysmax 
      real, pointer :: solnVar(:,:,:)

      solnVar(1:,1:,1:) => mainVar(:,:,:)

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
      nullify(solnVar) 

    end subroutine get_dt

    function compute_maxdivB(Bx, By) result(maxdivB)
      use sim_data, only: ilo, ihi, jlo, jhi, xTpts, yTpts, &
                            dx, dy
      implicit none
      real(8), intent(in), dimension(xTpts, yTpts) :: Bx, By
      real(8), dimension(xTpts, yTpts) :: divB
      real(8) :: maxdivB

      integer :: i, j
      
      divB = 0.0
      do i = ilo, ihi
        do j = jlo, jhi
           divB(i,j) = (Bx(i+1,j) - Bx(i-1,j)) / (2.0 * dx) &
                + (By(i,j+1) - By(i,j-1)) / (2.0 * dy)
           !if (divB > tol) then
           !  print *, "High error in divergence of B (> 1.0e-4)!"
           !  !print *, "Error at x, y = ", xvals(i), yvals(j)
           !  print *, "Error = ", divB
           !  stop
           !end if
        end do
      end do

      maxdivB = maxval(divB)
    end function compute_maxdivB

end module misc_module
