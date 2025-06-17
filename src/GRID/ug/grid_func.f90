module grid_func
#include "param.h"

  implicit none
  contains

  subroutine grid_init()
    use sim_data
    use misc_module, only: to_upper
    implicit none
  
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
   
    ! xrank and yrank 
    xrank = mod(myrank, xblk)
    yrank = floor(real(myrank)/xblk)
    
    ! default guard cells
    Gpts = 2  
    if (to_upper(trim(recon_method)) == "WENO3") Gpts = 2
    if (to_upper(trim(recon_method)) == "WENO5") Gpts = 3
 
    ! local interior points 
    lnx = nx / xblk
    lny = ny / yblk
  
    ! compute index of interior points 
    ilo = xrank * lnx + 1
    ihi = (xrank + 1) * lnx  
    jlo = yrank * lny + 1
    jhi = (yrank + 1) * lny 
 
    ! is the block at the boundary
    at_xlboundary = .false.
    at_xrboundary = .false.
    at_ylboundary = .false.
    at_yrboundary = .false.

    if (ilo == 1)  at_xlboundary = .true.
    if (ihi == nx) at_xrboundary = .true.
    if (jlo == 1)  at_ylboundary = .true.
    if (jhi == ny) at_yrboundary = .true.

    ! compute index of interior + guard cells
    iGlo = ilo - Gpts
    iGhi = ihi + Gpts
    jGlo = jlo - Gpts
    jGhi = jhi + Gpts
 
    ! total number of local points 
    lxTpts = lnx + 2 * Gpts
    lyTpts = lny + 2 * Gpts
  
    ! local min and max
    lxmin = xmin + (ilo - 1) * dx
    lxmax = xmin + ihi * dx
    lymin = ymin + (jlo - 1) * dx
    lymax = ymin + jhi * dx
  
    allocate(mainVar(NVAR_NUMBER, iGlo:iGhi, jGlo:jGhi))

    !set all the small values
    smalld = 1.0e-10
    smallp = 1.0e-10
    smalle = 1.0e-10
    small = 1.0e-100
 
    !print *, "myrank, xl, xr, yl, yr = ", myrank, at_xlboundary, at_xrboundary, &
    !                                              at_ylboundary, at_yrboundary 
    !print *, "myrank, ilo, ihi, jlo, jhi = ", myrank, ilo, ihi, jlo, jhi

  end subroutine grid_init

  subroutine grid_finalize()
    use sim_data, only: mainVar
    implicit none

    deallocate(mainVar)
  end subroutine

  subroutine get_coords(Dir,lo,hi,array)
    use sim_data, only: xmin, ymin, dx, dy
    use misc_module, only: to_upper 
    implicit none
    character(len=1), intent(in) :: Dir
    integer, intent(in) :: lo, hi
    real(8), dimension(lo:hi), intent(out) :: array
    integer :: ii
    real(8) :: delta, lmin
    
    array(lo:hi) = 0.0
 
    select case (to_upper(Dir))
    case ('X')
      lmin = xmin
      delta = dx
    case ('Y')
      lmin = ymin
      delta = dy
    case default
      print *, "Direction is X or Y only"
      stop
    end select

    do ii = lo, hi
      array(ii) = lmin + (2 * ii - 1) * delta/2.0
    end do 

  end subroutine  

end module grid_func
