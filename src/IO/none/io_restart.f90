module io_restart
#include "param.h"
  implicit none
contains
  subroutine restart_problem(fno, time)
    implicit none
    integer, intent(in) :: fno
    real, intent(out) :: time
    ! Do Nothing 
  end subroutine restart_problem
end module io_restart
