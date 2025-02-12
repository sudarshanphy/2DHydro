subroutine glm(bpsi, dt)
#include "header.h"

  use sim_data, only: ch, yTpts, xTpts
  implicit none
  real(8), dimension(xTpts, yTpts), intent(inout) :: bpsi
  real(8), intent(in) :: dt
  
  ! The ratio of the hyperbolic and parabolic effects
  !   c_r = ( c_p**2 / c_h ) ~ 0.18
  ! Ref: Dedner et al. (2002)  p.657, p.661

  real(8), parameter :: cr = 0.2d0
  real(8) :: exp_factor
  real(8), dimension(xTpts, yTpts) :: tmp_array
  
  !print *, "ch = ", ch
  exp_factor = exp( -0.5e0 * dt * ch / cr )
  
  tmp_array(:,:) = exp_factor * bpsi(:,:)

  bpsi(:,:) = tmp_array(:,:)

end subroutine glm
