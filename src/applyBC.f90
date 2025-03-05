module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applyBC(q, Dir, flip)

    use misc_module, only: to_upper
    use sim_data, only: xbctype, ybctype, &
                        iGlo, iGhi, jGlo, jGhi, Gpts
    implicit none
    real, pointer :: q(:, :)
    character(len=1), intent(in) :: Dir
    logical, optional, intent(in) :: flip
    real :: sig
    integer :: ii
    
    select case (to_upper(Dir))
    case ('X')
      if (to_upper(trim(xbctype)) == "PERIODIC") then
        do ii  = iGlo, ilo - 1
          ! lower face
          q(ii, :) = q(iGhi-2*Gpts+ii, :)
          ! upper face
          q(iGhi-Gpts+ii, :) = q(Gpts + ii, :)
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
    use sim_data, only: mainVar, iGlo, iGhi, &
                        jGlo, jGhi, at_xboundary, &
                        at_yboundary
    implicit none
    real, pointer :: q(:,:)
    integer :: n

    if (at_xboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(:,:,n)
          if (n == VELX_VAR) then
            call applyBC(q, "x", .true.)
          else if (n == VELY_VAR) then
            call applyBC(q,"x")
          else
            call applyBC(q, "x")
         end if
         nullify(q) 
       end do
    endif

    if (at_yboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(:,:,n)
          if (n == VELX_VAR) then
            call applyBC(q,"y")
          else if (n == VELY_VAR) then
            call applyBC(q, "y", .true.)
          else
            call applyBC(q, "y")
         end if
         nullify(q) 
       end do
    endif
    
  end subroutine applyBC_all
end module applyBC_module

