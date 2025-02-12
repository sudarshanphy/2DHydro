module io_module
#include "header.h"

  implicit none
contains
  
  subroutine write_output(t, step, x, y, &
                          d, u, v, p, e, &
#ifdef MHD 
                          bx, by, psi,   &
#endif
                          &outputno, restart_init)
   use sim_data, only: nx, ny, gamma, basenm, dx, dy
   implicit none
  real, intent(in) :: t
  integer, intent(in) :: step, outputno
  real, dimension(nx), intent(in) :: x
  real, dimension(ny), intent(in) :: y
  real, dimension(nx, ny), intent(in) :: d, u, v, p, e
#ifdef MHD
  real, dimension(nx, ny), intent(in) :: bx, by, psi
#endif
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
  write(ionum, *) "# xcenter  ycenter  dens  velx  vely  pres  ener  bmfx  bmfy  bpsi"
#else
  write(ionum, *) "# xcenter  ycenter  dens  velx  vely  pres  ener"
#endif
  do i = 1, nx
     do j = 1, ny
#ifdef MHD
        write(ionum, 12) x(i), y(j), d(i,j), u(i,j), v(i,j), p(i,j), e(i,j), &
                                     &  bx(i,j), by(i,j), psi(i,j)
#else
        write(ionum, 12) x(i), y(j), d(i,j), u(i,j), v(i,j), p(i,j), e(i,j)
#endif
     end do
     write(ionum, *) "####"
  end do
  close(ionum)
  end subroutine write_output
end module io_module
