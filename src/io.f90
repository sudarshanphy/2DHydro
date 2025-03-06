module io_module
#include "param.h"

  implicit none
contains
  
  subroutine write_output(t, step, outputno, restart_init)
   use sim_data, only: nx, ny, gamma, basenm, dx, dy, &
                       mainVar, ilo, ihi, jlo, jhi, &
                       Gpts, myrank, xblk, yblk, lnx, lny

   use grid_func, only: get_coords
   implicit none
   real, intent(in) :: t
   integer, intent(in) :: step, outputno
   real, pointer :: solnVar(:,:,:)
   logical, optional, intent(in) :: restart_init
   integer :: ionum, i, j
   character(len=256) :: fname
   character(len=10) :: int_to_str, rank_to_str

   real, allocatable:: x(:), y(:)

  allocate(x(ilo:ihi))
  allocate(y(jlo:jhi))

  call get_coords('x',ilo,ihi,x)
  call get_coords('y',jlo,jhi,y)

  solnVar(ilo:,jlo:,1:) => mainVar(ilo:ihi, jlo:jhi, 1:)
  ionum = 100
  
12   format (1x, 50(es25.18, :, 1x))
  
  write(int_to_str, "(I4.4)") outputno
  write(rank_to_str, "(I4.4)") myrank

  if (present(restart_init)) then
    fname = "./output/"//trim(adjustl(basenm))//"_reinit_"//trim(adjustl(rank_to_str))//&
      "_"//trim(adjustl(int_to_str))//".dat"
  else
    fname = "./output/"//trim(adjustl(basenm))//"_"//trim(adjustl(rank_to_str))//&
      "_"//trim(adjustl(int_to_str))//".dat"
  endif

  open(unit=ionum, file=fname, status="replace")
  write(ionum, *) "#time = ", t, ", step = ", step
  write(ionum, *) "#nx = ", nx , ", dx = ", dx
  write(ionum, *) "#ny = ", ny , ", dy = ", dy
  write(ionum, *) "#lnx = ", lnx, ", xblk = ", xblk
  write(ionum, *) "#lny = ", lny, ", yblk = ", yblk
  write(ionum, *) "#gamma = ", gamma
#ifdef MHD
  write(ionum, *) "# xcenter ycenter dens velx vely velz pres bmfx bmfy bmfz bpsi ener"
#else
  write(ionum, *) "# xcenter ycenter dens velx vely velz pres ener"
#endif
  do i = ilo, ihi
     do j = jlo, jhi
        write(ionum, 12) x(i), y(j), solnVar(i,j,:)
     end do
     write(ionum, *) "####"
  end do
  close(ionum)
  deallocate(x,y)
  nullify(solnVar)
  end subroutine write_output
end module io_module
