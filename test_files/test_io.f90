program io_test
  use io_module, only: write_output

  implicit none
  integer :: nx, ny, stepno, outputno
  real :: time, gamma
  real, allocatable, dimension(:) :: x, y
  real, allocatable, dimension(:, :) :: d, u, v, p, e
  character(len=256) :: basenm
  real :: dx, dy, xmin, xmax, ymin, ymax
  integer :: i, j, ilo, ihi, jlo, jhi

  nx = 9; ny = 7; stepno = 10; outputno = 4
  time = 0.05; gamma = 1.4
  xmin = 0.0; ymin = 0.0
  xmax = 1.0; ymax = 1.0
  dx = (xmax - xmin)/(nx-4)
  dy = (ymax - ymin)/(ny-4)
  ilo = 3; ihi = nx - 2; jlo = 3; jhi = ny - 2

  allocate(x(nx), y(ny))
  allocate(d(nx, ny), u(nx, ny), v(nx, ny), p(nx, ny), e(nx, ny)) 

  do j = 1, ny
     do i = 1, nx
        x(i) = xmin - 2.0d0 * dx +  (2 * i  - 1)* dx / 2
        y(j) = ymin - 2.0d0 * dy + (2 * j - 1) * dy / 2
        d(i, j) = 5.0 + i + j
        u(i, j) = 6.0 + i + j
        v(i, j) = 7.0 + i + j
        p(i, j) = 8.0 + i + j
        e(i, j) = 9.0 + i + j
     end do
  end do
  basenm = "io_test"
  print *, x
  print *, y
  call write_output(time, stepno, nx-4, ny-4, x(ilo:ihi), y(jlo:jhi), &
    d(ilo:ihi, jlo:jhi), u(ilo:ihi, jlo:jhi), v(ilo:ihi, jlo:jhi), &
    p(ilo:ihi, jlo:jhi), e(ilo:ihi, jlo:jhi), gamma, basenm, outputno)

  deallocate(x, y, d, u, v, p, e)
end program io_test
