module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applydefBC(var, Dir, Face, bc_int)

    use sim_data, only: xlbc_int, ylbc_int, &
                        xrbc_int, yrbc_int, &
                        Gpts, ihi, ilo, jhi, jlo, &
                        iGlo, jGlo, iGhi, jGhi, myrank, &
                        mainVar
    implicit none
    real, pointer :: q(:,:)
    integer, intent(in) :: var, Dir, Face, bc_int
    real :: sig
    integer :: ii

    if (bc_int == CUSTOM) return

    q(iGlo:, jGlo:) => mainVar(var,iGlo:iGhi,jGlo:jGhi)

    if (Face == LEFT .and. Dir == IAXIS) then
      if (bc_int == FLOW) then
        do ii = 1, Gpts
          q(ilo-ii, :) = q(ilo, :)
        end do
      else if (bc_int == REFLECT) then
        sig = 1.000
        if (var == VELX_VAR) sig = -1.000
        do ii = 1, Gpts
          q(ilo-ii, :) = sig * q(ilo+ii-1, :)
        end do
      end if
    else if (Face == LEFT .and. Dir == JAXIS) then
      if (bc_int == FLOW) then
        do ii = 1, Gpts
          q(:, jlo-ii) = q(:, jlo)
        end do
      else if (bc_int == REFLECT) then
        sig = 1.000
        if (var == VELY_VAR) sig = -1.000
        do ii = 1, Gpts
          q(:, jlo-ii) = sig * q(:, jlo+ii-1)
        end do
      end if
    else if (Face == RIGHT .and. Dir == IAXIS) then
      if (bc_int == FLOW) then
        do ii = 1, Gpts
          q(ihi+ii, :) = q(ihi, :)
        end do
      else if (bc_int == REFLECT) then
        sig = 1.000
        if (var == VELX_VAR) sig = -1.000
        do ii = 1, Gpts
          q(ihi+ii, :) = sig * q(ihi-ii+1, :)
        end do
      end if
    else if (Face == RIGHT .and. Dir == JAXIS) then
      if (bc_int == FLOW) then
        do ii = 1, Gpts
          q(:, jhi+ii) = q(:, jhi)
        end do
      else if (bc_int == REFLECT) then
        sig = 1.000
        if (var == VELY_VAR) then
                sig = -1.000
        end if
        do ii = 1, Gpts
          q(:, jhi+ii) = sig * q(:, jhi-ii+1)
        end do
      end if
    else
      print *, "Wrong!! This code only solves in 2D"
      stop
    end if

    nullify(q)
  end subroutine applydefBC

  subroutine applyBC_all()
    use sim_data, only: at_xlboundary, &
                        at_ylboundary, at_xrboundary, &
                        at_yrboundary, xlbc_int, xrbc_int, &
                        ylbc_int, yrbc_int 
    use customBC_module, only: applycustomBC
    implicit none
    integer :: n
    
    ! Periodic BC is already applied in the guardcell_fill routine
    ! applyBC_all should always be called after guardcell_fill routine

    if (at_xlboundary) then
       if (xlbc_int == CUSTOM) then
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applycustomBC(n,IAXIS,LEFT)
          end do

       else 
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applydefBC(n,IAXIS,LEFT,xlbc_int)
          end do
       end if
    endif

    if (at_xrboundary) then
       if (xrbc_int == CUSTOM) then
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applycustomBC(n,IAXIS,RIGHT)
          end do
       
       else 
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applydefBC(n,IAXIS,RIGHT,xrbc_int)
          end do
       end if
    endif

    if (at_ylboundary) then
       if (ylbc_int == CUSTOM) then
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applycustomBC(n,JAXIS,LEFT)
          end do
       
       else 
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applydefBC(n,JAXIS,LEFT,ylbc_int)
          end do
       end if
    endif
    
    if (at_yrboundary) then 
       if (yrbc_int == CUSTOM) then
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applycustomBC(n,JAXIS,RIGHT)
          end do
       
       else
          do n=NVAR_BEGIN, NVAR_NUMBER
             call applydefBC(n,JAXIS,RIGHT,yrbc_int)
          end do
       end if
    endif

  end subroutine applyBC_all
end module applyBC_module

