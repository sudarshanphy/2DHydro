module io_module
#include "param.h"

  implicit none
contains
  
  subroutine write_output(t, step, outputno, restart_init)
      implicit none
      real, intent(in) :: t
      integer, intent(in) :: step, outputno
      real, pointer :: solnVar(:,:,:)
      logical, optional, intent(in) :: restart_init

      ! Do nothing! -> No outputfiles
  end subroutine write_output
end module io_module
