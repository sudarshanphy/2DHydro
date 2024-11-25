module misc_module
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

    function get_dt(dens, velx, vely, pres) result(dt)
      use sim_data, only: gamma, xTpts, yTpts, dx, dy, cfl
      implicit none
      real, dimension(xTpts, yTpts) :: dens, velx, vely, pres
      real, dimension(xTpts, yTpts) :: cs, xsmax, ysmax
      !real :: cs, xsmax_init, ysmax_init, xsmax, ysmax
      real :: dt
      integer :: i, j
      real :: max_xsmax, max_ysmax 
      !xsmax_init = 1.0e-9
      !ysmax_init = 1.0e-9 
      !do j = 1, yTpts
      !  do i = 1, xTpts
      !     cs = sqrt(gamma * pres(i,j) / dens(i,j))
      !     xsmax = max(cs + abs(velx(i,j)), xsmax_init)
      !     ysmax = max(cs + abs(vely(i,j)), ysmax_init)
      !     xsmax_init = xsmax
      !     ysmax_init = ysmax
      !  end do
      !end do
      !dt = min(cfl * dx/xsmax, cfl * dy/ysmax) 
      
      do j = 1, yTpts
        do i = 1, xTpts
           cs(i,j) = sqrt(gamma * pres(i,j) / dens(i,j))
           xsmax(i,j) = cs(i,j) + abs(velx(i,j))
           ysmax(i,j) = cs(i,j) + abs(vely(i,j))
        end do
      end do
       
      max_xsmax = maxval(xsmax); max_ysmax = maxval(ysmax)

      !print *, "max xsmax, ysmax = ", max_xsmax, max_ysmax 
      dt = min(cfl * dx/max_xsmax, cfl * dy/max_ysmax) 

    end function get_dt

end module misc_module
