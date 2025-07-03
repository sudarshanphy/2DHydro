module customBC_module
#include "param.h"
  implicit none
contains

  subroutine applycustomBC(var, dir, face)
                        
    use sim_data, only: Gpts, ihi, ilo, jhi, jlo, &
                        jGlo, jGhi, iGlo, iGhi, PI, &
                        myrank, mainVar
    use grid_func, only: get_coords
    implicit none
    real, pointer :: q(:, :)
    integer, intent(in) :: var, Dir, Face
    real :: sig
    integer :: ii, j, i
    real, allocatable, dimension(:) :: x,y

    !allocate(x(iGlo:iGhi), y(jGlo:jGhi))
    ! 
    !q(iGlo:, jGlo:) => mainVar(var,:,:)

    !call get_coords(IAXIS,iGlo,iGhi,x)
    !call get_coords(JAXIS,jGlo,jGhi,y)

    !if ((Dir == IAXIS) .and. (Face == LEFT)) then
    !    ! Apply boundary at XL
    !end if
    ! 
    !nullify(q) 
    !deallocate(x,y)
  end subroutine applycustomBC

end module customBC_module
