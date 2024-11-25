module io_module
  implicit none
contains
  
  subroutine write_output(t, step, x, y, &
                          d, u, v, p, e, outputno)
   use sim_data, only: nx, ny, gamma, basenm, dx, dy
   implicit none
  real, intent(in) :: t
  integer, intent(in) :: step, outputno
  real, dimension(nx), intent(in) :: x
  real, dimension(ny), intent(in) :: y
  real, dimension(nx, ny), intent(in) :: d, u, v, p, e
  integer :: ionum, i, j
  character(len=256) :: fname
  character(len=10) :: int_to_str

  ionum = 100

12   format (1x, 50(es25.18, :, 1x))
  
  write(int_to_str, "(I4.4)") outputno
  fname = trim(adjustl(basenm))//"_"//trim(adjustl(int_to_str))//".dat"
  open(unit=ionum, file=fname, status="replace")
  write(ionum, *) "#time = ", t, ", step = ", step
  write(ionum, *) "#nx = ", nx , ", dx = ", dx
  write(ionum, *) "#ny = ", ny , ", dy = ", dy
  write(ionum, *) "#gamma = ", gamma
  write(ionum, *) "# xcenter  ycenter  dens  velx  vely  pres  ener"
  do i = 1, nx
     do j = 1, ny
        write(ionum, 12) x(i), y(j), d(i,j), u(i,j), v(i,j), p(i,j), e(i,j)
     end do
     write(ionum, *) "####"
  end do
  close(ionum)
  end subroutine write_output
end module io_module
