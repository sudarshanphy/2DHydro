module io_module
#include "param.h"

  implicit none
contains
  
  subroutine write_output(t, step, x, y, solnVar,&
                           outputno, restart_init)
   use sim_data, only: nx, ny, gamma, basenm, dx, dy
   implicit none
  real, intent(in) :: t
  integer, intent(in) :: step, outputno
  real, dimension(nx), intent(in) :: x
  real, dimension(ny), intent(in) :: y
  real, dimension(nx, ny, NVAR_NUMBER), intent(in) :: solnVar
  logical, optional, intent(in) :: restart_init
  integer :: ionum, i, j
  character(len=256) :: fname
  character(len=10) :: int_to_str

  ionum = 100

12   format (1x, 50(es25.18, :, 1x))
  
  write(int_to_str, "(I4.4)") outputno
  if (present(restart_init)) then
    fname = "./output/"//trim(adjustl(basenm))//"_reinit_"//trim(adjustl(int_to_str))//".dat"
  else
    fname = "./output/"//trim(adjustl(basenm))//"_"//trim(adjustl(int_to_str))//".dat"
  endif

  open(unit=ionum, file=fname, status="replace")
  write(ionum, *) "#time = ", t, ", step = ", step
  write(ionum, *) "#nx = ", nx , ", dx = ", dx
  write(ionum, *) "#ny = ", ny , ", dy = ", dy
  write(ionum, *) "#gamma = ", gamma
#ifdef MHD
  write(ionum, *) "# xcenter ycenter dens velx vely velz pres bmfx bmfy bmfz bpsi ener"
#else
  write(ionum, *) "# xcenter ycenter dens velx vely velz pres ener"
#endif
  do i = 1, nx
     do j = 1, ny
        write(ionum, 12) x(i), y(j), solnVar(i,j,:)
     end do
     write(ionum, *) "####"
  end do
  close(ionum)
  end subroutine write_output
end module io_module
