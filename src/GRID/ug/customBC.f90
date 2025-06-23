module customBC
#include "param.h"
  implicit none
contains

  subroutine applycustomBC(var, q, dir, face, notapplied)
    use sim_data, only: xlbc_int, ylbc_int, &
                        xrbc_int, yrbc_int, &
                        Gpts, ihi, ilo, jhi, jlo
    implicit none
    real, pointer :: q(:, :)
    integer, intent(in) :: var, Dir, Face
    logical, intent(inout) :: notapplied
    real :: sig
    integer :: ii
    
    select case (Face)
    ! BC on the left face
    case (LEFT)

      select case (Dir)
      !BC in X-direction
      case (IAXIS)
        !if (to_upper(trim(xlbctype)) == "CUSTOM") then
           ! Apply boundary at XL
      case (JAXIS)
           ! Apply boundary at YL
      case default

      end select

    case (RIGHT)

      select case (Dir)
      !BC in X-direction
      case (IAXIS)
        !if (to_upper(trim(xlbctype)) == "CUSTOM") then
           ! Apply boundary at XR
      case (JAXIS)
           ! Apply boundary at YR
      case default

      end select
    end select

  end subroutine applycustomBC

end module customBC
