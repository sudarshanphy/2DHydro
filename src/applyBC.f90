module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applyBC(q, Dir, Face, flip)

    use misc_module, only: to_upper
    use sim_data, only: xlbctype, ylbctype, &
                        xrbctype, yrbctype, &
                        Gpts, ihi, ilo, jhi, jlo
    implicit none
    real, pointer :: q(:, :)
    character(len=1), intent(in) :: Dir, Face
    logical, optional, intent(in) :: flip
    real :: sig
    integer :: ii
    
    select case (to_upper(Face))
    ! BC on the left face
    case ('L')

      select case (to_upper(Dir))
      !BC in X-direction
      case ('X')
        if (to_upper(trim(xlbctype)) == "PERIODIC") then
          do ii  = 1, Gpts
            ! lower face
            !q(ii, :) = q(iGhi-2*Gpts+ii, :)
          end do

        else if (to_upper(trim(xlbctype)) == "FLOW") then
          do ii  = 1, Gpts
            ! lower face
            q(ilo-ii, :) = q(ilo, :)
          end do

        else if (to_upper(trim(xlbctype)) == "REFLECT") then
          sig = 1.000
          if (present(flip)) sig = -1.000
          do ii = 1, Gpts
            ! lower face
            q(ilo-ii, :) = sig * q(ilo+ii-1, :)
          end do
        end if

      ! BC in Y-direction  
      case ('Y')
        if (to_upper(trim(ylbctype)) == "PERIODIC") then
          do ii  = 1, Gpts
            ! lower face
            !q(:, ii) = q(:, yTpts-2*Gpts+ii)
          end do

        else if (to_upper(trim(ylbctype)) == "FLOW") then
          do ii  = 1, Gpts
            ! lower face
            q(:, jlo-ii) = q(:, jlo)
          end do

        else if (to_upper(trim(ylbctype)) == "REFLECT") then
          sig = 1.000
          if (present(flip)) sig = -1.000
          do ii = 1, Gpts
            ! lower face
            q(:, jlo-ii) = sig * q(:, jlo+ii-1)
          end do
        end if

      case default
          print *, "Wrong!! This code only solves in 2D"
          stop
      end select
    
    ! BC on the right face
    case ('R')
      select case (to_upper(Dir))
      ! BC in X-direction
      case ('X')
        if (to_upper(trim(xrbctype)) == "PERIODIC") then
          do ii  = 1, Gpts
            ! upper face
            !q(iGhi-Gpts+ii, :) = q(Gpts + ii, :)
          end do

        else if (to_upper(trim(xrbctype)) == "FLOW") then
          do ii  = 1, Gpts
            ! upper face
            q(ihi+ii, :) = q(ihi, :)
          end do

        else if (to_upper(trim(xrbctype)) == "REFLECT") then
          sig = 1.000
          if (present(flip)) sig = -1.000
          do ii = 1, Gpts
            ! upper face
            q(ihi+ii, :) = sig * q(ihi-ii+1, :)
          end do
        end if

      ! BC in Y-direction
      case ('Y')
        if (to_upper(trim(yrbctype)) == "PERIODIC") then
          do ii  = 1, Gpts
            ! upper face
            !q(:, yTpts-Gpts+ii) = q(:, Gpts + ii)
          end do

        else if (to_upper(trim(yrbctype)) == "FLOW") then
          do ii  = 1, Gpts
            ! upper face
            q(:, jhi+ii) = q(:, jhi)
          end do

        else if (to_upper(trim(yrbctype)) == "REFLECT") then
          sig = 1.000
          if (present(flip)) sig = -1.000
          do ii = 1, Gpts
            ! upper face
            q(:, jhi+ii) = sig * q(:, jhi-ii+1)
          end do
        end if 

      case default
          print *, "Wrong!! This code only solves in 2D"
          stop
      end select

    end select  
  end subroutine applyBC
  
  subroutine applyBC_all()
    use sim_data, only: mainVar, iGlo, iGhi, &
                        jGlo, jGhi, at_xlboundary, &
                        at_ylboundary, at_xrboundary, &
                        at_yrboundary
    implicit none
    real, pointer :: q(:,:)
    integer :: n

    if (at_xlboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          if (n == VELX_VAR) then
            call applyBC(q,"x","L", .true.)
          else
            call applyBC(q,"x","L")
         end if
         nullify(q) 
       end do
    endif

    if (at_xrboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          if (n == VELX_VAR) then
            call applyBC(q,"x","R", .true.)
          else
            call applyBC(q,"x","R")
         end if
         nullify(q) 
       end do
    endif

    if (at_ylboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          if (n == VELY_VAR) then
            call applyBC(q,"y","L", .true.)
          else
            call applyBC(q,"y","L")
         end if
         nullify(q) 
       end do
    endif
    
    if (at_yrboundary) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          if (n == VELY_VAR) then
            call applyBC(q,"y","R", .true.)
          else
            call applyBC(q,"y","R")
         end if
         nullify(q) 
       end do
    endif

  end subroutine applyBC_all
end module applyBC_module

