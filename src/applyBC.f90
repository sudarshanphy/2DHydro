module applyBC_module
#include "header.h"
#include "param.h"
  implicit none
contains

  subroutine applyBC(q, Dir, flip)

    use misc_module, only: to_upper
    use sim_data, only: xbctype, ybctype, &
                        xTpts, yTpts, Gpts
    implicit none
    real(8), dimension(xTpts, yTpts), intent(inout) :: q
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
  
  subroutine applyBC_all(solnVar)
    use sim_data, only: xTpts, yTpts
    implicit none
    real, dimension(xTpts, yTpts, NVAR_NUMBER), intent(inout) :: solnVar
    integer :: n

    do n=NVAR_BEGIN, NVAR_NUMBER
       if (n == VELX_VAR) then
         call applyBC(solnVar(:,:,n), "x", .true.)
         call applyBC(solnVar(:,:,n),"y")
       else if (n == VELY_VAR) then
         call applyBC(solnVar(:,:,n), "y", .true.)
         call applyBC(solnVar(:,:,n),"x")
       else
         call applyBC(solnVar(:,:,n), "x")
         call applyBC(solnVar(:,:,n), "y")
      end if
    end do
     
  end subroutine applyBC_all
end module applyBC_module

