module customBC
#include "param.h"
  implicit none
contains

  subroutine applycustomBC(var, q, dir, face)
    use sim_data, only: xlbc_int, ylbc_int, &
                        xrbc_int, yrbc_int, &
                        Gpts, ihi, ilo, jhi, jlo, &
                        jGlo, jGhi, iGlo, iGhi, PI
    use grid_func, only: get_coords
    implicit none
    real, pointer :: q(:, :)
    integer, intent(in) :: var, Dir, Face
    real :: sig
    integer :: ii, j, i
    real, allocatable, dimension(:) :: x,y

    allocate(x(iGlo:iGhi), y(jGlo:jGhi))
    
    call get_coords(IAXIS,iGlo,iGhi,x)
    call get_coords(JAXIS,jGlo,jGhi,y)

    select case (Face)
    ! BC on the left face
    case (LEFT)

      select case (Dir)
      !BC in X-direction
      case (IAXIS)
           ! Apply boundary at XL
          sig = 1.000
          if (var == VELX_VAR) sig = -1.000
           ! reflecting boundary
          do ii = 1, Gpts
            ! lower face
            q(ilo-ii, :) = sig * q(ilo+ii-1, :)
          end do
          ! inside a certain region
          do j = jGlo, jGhi
            do i = iGlo, iGhi
              if (abs(y(j)) <= 0.01) then
                if (x(i) <= 0.25) then
                   if (var==VELX_VAR) q(i,j) = 0.5
                   if (var==DENS_VAR) q(i,j) = 2.0 
                end if
              end if
            end do
          end do

      case (JAXIS)
           ! Apply boundary at YL
      case default

      end select

    case (RIGHT)

      select case (Dir)
      !BC in X-direction
      case (IAXIS)
           ! Apply boundary at XR
      case (JAXIS)
           ! Apply boundary at YR
      case default

      end select
    end select
    
    deallocate(y)
  end subroutine applycustomBC

end module customBC
