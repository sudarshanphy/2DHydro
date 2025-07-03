subroutine glm(dt)
#include "param.h"

  use sim_data, only: ch, mainVar, iGlo, jGlo, &
                                   iGhi, jGhi
  implicit none
  real, pointer :: bpsi(:,:)
  real(8), intent(in) :: dt
  
  ! The ratio of the hyperbolic and parabolic effects
  !   c_r = ( c_p**2 / c_h ) ~ 0.18
  ! Ref: Dedner et al. (2002)  p.657, p.661

  real(8), parameter :: cr = 0.2d0
  real(8) :: exp_factor
  real(8), dimension(iGlo:iGhi, jGlo:jGhi) :: tmp_array

#ifdef MHD  
  bpsi(iGlo:,jGlo:) => mainVar(BPSI_VAR,:,:)
  !print *, "ch = ", ch
  exp_factor = exp( -0.5e0 * dt * ch / cr )
  
  tmp_array(:,:) = exp_factor * bpsi(:,:)

  bpsi(:,:) = tmp_array(:,:)
  
  nullify(bpsi)
#endif
end subroutine glm
