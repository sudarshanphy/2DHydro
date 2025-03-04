module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applyBC(q, Dir, flip)

    use misc_module, only: to_upper
    use sim_data, only: xbctype, ybctype, &
                        xTpts, yTpts, Gpts
    implicit none
    real, pointer :: q(:, :)
    character(len=1), intent(in) :: Dir
    logical, optional, intent(in) :: flip
    real :: sig
    integer :: ii
    
    select case (to_upper(Dir))
    case ('X')
      if (to_upper(trim(xbctype)) == "PERIODIC") then
        do ii  = 1, Gpts
          ! lower face
          q(ii, :) = q(xTpts-2*Gpts+ii, :)
          ! upper face
          q(xTpts-Gpts+ii, :) = q(Gpts + ii, :)
        end do
      else if (to_upper(trim(xbctype)) == "FLOW") then
        do ii  = 1, Gpts
          ! lower face
          q(ii, :) = q(Gpts+1, :)
          ! upper face
          q(xTpts-Gpts+ii, :) = q(xTpts-Gpts, :)
        end do
      else if (to_upper(trim(xbctype)) == "REFLECT") then
        sig = 1.000
        if (present(flip)) sig = -1.000
        do ii = 1, Gpts
          ! lower face
          q(ii, :) = sig * q(2*Gpts+1-ii, :)
          ! upper face
          q(xTpts-Gpts+ii, :) = sig * q(xTpts-Gpts+1-ii, :)
        end do
      end if 
    case ('Y')
      if (to_upper(trim(ybctype)) == "PERIODIC") then
        do ii  = 1, Gpts
          ! lower face
          q(:, ii) = q(:, yTpts-2*Gpts+ii)
          ! upper face
          q(:, yTpts-Gpts+ii) = q(:, Gpts + ii)
        end do
      else if (to_upper(trim(ybctype)) == "FLOW") then
        do ii  = 1, Gpts
          ! lower face
          q(:, ii) = q(:, Gpts+1)
          ! upper face
          q(:, yTpts-Gpts+ii) = q(:, ytpts-Gpts)
        end do
      else if (to_upper(trim(ybctype)) == "REFLECT") then
        sig = 1.000
        if (present(flip)) sig = -1.000
        do ii = 1, Gpts
          ! lower face
          q(:, ii) = sig * q(:, 2*Gpts+1-ii)
          ! upper face
          q(:, yTpts-Gpts+ii) = sig * q(:, yTpts-Gpts+1-ii)
        end do
      end if 
    case default
        print *, "Wrong!! This code only solves in 2D"
        stop
    end select
  end subroutine applyBC
  
  subroutine applyBC_all()
    use sim_data, only: xTpts, yTpts, mainVar
    implicit none
    real, pointer :: q(:,:)
    integer :: n


    do n=NVAR_BEGIN, NVAR_NUMBER
       q(1:,1:) => mainVar(1:,1:,n)
       if (n == VELX_VAR) then
         call applyBC(q, "x", .true.)
         call applyBC(q,"y")
       else if (n == VELY_VAR) then
         call applyBC(q, "y", .true.)
         call applyBC(q,"x")
       else
         call applyBC(q, "x")
         call applyBC(q, "y")
      end if
      nullify(q) 
    end do
    
  end subroutine applyBC_all
end module applyBC_module

