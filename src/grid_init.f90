subroutine grid_init(xIcenter, yIcenter)

  use sim_data
  use misc_module, only: to_upper
  implicit none
  real, intent(inout), dimension(nx) :: xIcenter
  real, intent(inout), dimension(ny) :: yIcenter
  integer :: i, j

  dx = (xmax - xmin) / nx
  dy = (ymax - ymin) / ny
  
  if (to_upper(trim(recon_method)) == "WENO3") Gpts = 2
  if (to_upper(trim(recon_method)) == "WENO5") Gpts = 3
  
  xTpts = nx + 2 * Gpts
  yTpts = ny + 2 * Gpts

  ilo = Gpts + 1; jlo = Gpts + 1
  ihi = nx + Gpts; jhi = ny + Gpts

  do i =1, nx
    xIcenter(i) = xmin + (2 * i - 1) * dx / 2.0d0
  end do 
  do j =1, ny
    yIcenter(j) = ymin + (2 * j - 1) * dy / 2.0d0
  end do 

end subroutine grid_init
