module resistiveflux
#include "param.h"
  implicit none
  contains

    function eta(xx,yy) result(val)
      implicit none
      real(8), intent(in) :: xx, yy
      real(8) :: val
      real(8) :: eta0, eta01
      ! Resistivity
      real(8), parameter :: Rm1 = 60.0e0
      real(8), parameter :: Rm0 = 1000.0e0

      eta0 = 1.0e0/Rm0
      eta01 = (Rm0 - Rm1)/ (Rm0 * Rm1)
      val = eta0 + eta01 * ( cosh(min( sqrt( xx**2+yy**2 ), 25.0e0 )) )**(-2)

    end function eta


    subroutine get_resistive_dt(dt)
      use sim_data, only: cfl, dx, dy
      implicit none
      real(8), intent(inout) :: dt
      real(8), parameter :: Rm1 = 60.0e0

      dt = min( dt, min(0.5d0*cfl*Rm1*(dx**2), 0.5d0*cfl*Rm1*(dy**2)) )

    end subroutine get_resistive_dt

#ifdef MHD

    subroutine resistive_Flux(U, F,i,j,Dir)
      use sim_data, only: ilo, ihi, jlo, jhi, xTpts, yTpts, &
                          dx, dy, xmin, ymin, Gpts
      use misc_module, only: to_upper
      implicit none
      integer, intent(in) :: i, j
      real(8), dimension(xTpts, yTpts), intent(in):: U
      real(8), dimension(xTpts, yTpts), intent(inout):: F
      character(len=1), intent(in) :: Dir
      real(8) :: xf, yf  ! x and y coords at faces
      real(8) :: etaf, delta, Jt1, Jt2
      real(8), dimension(2) :: n
      integer :: bn, bt1, bt2

      select case (to_upper(Dir))
      case ('X')
          n = (/1.00, 0.00/)
          delta = 1.0 / dx
          bn  = BMFX_VAR
          ! tangential magnetic fields
          bt1 = BMFY_VAR
          bt2 = BMFZ_VAR
      case ('Y')
          n = (/0.00, 1.00/)
          delta = 1.0 / dy
          bn  = BMFY_VAR
          bt1 = BMFX_VAR
          bt2 = BMFZ_VAR
      case default
          print *, "Wrong!! Direction should be X or Y!"
          stop
      end select

      xf = xmin + (i - (Gpts + 1)) * dx + n(2) * dx/2
      yf = ymin + (j - (Gpts + 1)) * dy + n(1) * dy/2

      etaf = eta(xf, yf) * delta

      !Jt1 = -etaf * ()

    end subroutine resistive_Flux

#endif
end module resistiveflux
