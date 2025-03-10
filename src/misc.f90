module misc_module
#include "param.h"
  use mpi
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
      use sim_data, only: gamma, dx, dy, cfl, &
                          mainVar, iGlo, iGhi, &
                          jGlo, jGhi, comm, ierr
      use mpi
#ifdef MHD
      use sim_data, only: ch
#endif
      implicit none
      real, dimension(iGlo:iGhi, jGlo:jGhi) :: xsmax, ysmax
      real :: cs, xcmax, ycmax   !sound speed
#ifdef MHD  
      real :: cax, cay, cfx, cfy, B2, cB2 !alfven wave and magnetosonic wave
#endif
      real, intent(out) :: dt
      real :: localdt, localch
      integer :: i, j
      real :: max_xsmax, max_ysmax 
      real, pointer :: solnVar(:,:,:)

      solnVar(1:,iGlo:,jGlo:) => mainVar(:,:,:)

      do j = jGlo, jGhi
        do i = iGlo, iGhi
           cs = sqrt(gamma * solnVar(PRES_VAR,i,j) / solnVar(DENS_VAR,i,j))

           if (cs < 0.0) then
             print *, "Imaginary sound speed!"
             print *, "At i,j, with gamma, pres, dens = ", i, j, gamma, &
                           solnVar(PRES_VAR,i,j), solnVar(DENS_VAR,i,j)
           endif
           xcmax = cs
           ycmax = cs
#ifdef MHD
           cax = solnVar(BMFX_VAR,i,j) / sqrt(solnVar(DENS_VAR,i,j))
           cay = solnVar(BMFY_VAR,i,j) / sqrt(solnVar(DENS_VAR,i,j))
           B2 =  sum(solnVar(BMFX_VAR:BMFZ_VAR,i,j)*solnVar(BMFX_VAR:BMFZ_VAR,i,j))
           cB2 = B2 / solnVar(DENS_VAR,i,j)
           cfx = sqrt(0.5 * ((cs + cB2) + sqrt((cs + cB2)**2 - 4.0 * cs * cax * cax)))
           cfy = sqrt(0.5 * ((cs + cB2) + sqrt((cs + cB2)**2 - 4.0 * cs * cay * cay)))
           xcmax = cfx
           xcmax = cfy
#endif
           xsmax(i,j) = xcmax + abs(solnVar(VELX_VAR,i,j))
           ysmax(i,j) = ycmax + abs(solnVar(VELY_VAR,i,j))
        end do
      end do
       
      max_xsmax = maxval(xsmax); max_ysmax = maxval(ysmax)
      
#ifdef MHD
      !speed for divergence cleaning
      localch = max(max_xsmax, max_ysmax)
      ! get max ch from all the cores
      call MPI_ALLREDUCE(localch, ch, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)
#endif

      !print *, "max xsmax, ysmax = ", max_xsmax, max_ysmax 
      localdt = min(cfl * dx/max_xsmax, cfl * dy/max_ysmax)
     
      !print *, "local dt = ", localdt 
      ! get min dt from all the cores
      call MPI_ALLREDUCE(localdt, dt, 1, MPI_DOUBLE, MPI_MIN, comm, ierr)
      !print *, "global dt = ", dt
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
