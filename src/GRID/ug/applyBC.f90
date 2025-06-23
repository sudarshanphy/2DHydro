module applyBC_module
#include "param.h"
  implicit none
contains

  subroutine applyBC(var, q, Dir, Face)

    use sim_data, only: xlbc_int, ylbc_int, &
                        xrbc_int, yrbc_int, &
                        Gpts, ihi, ilo, jhi, jlo
    implicit none
    real, pointer :: q(:, :)
    integer, intent(in) :: var, Dir, Face
    real :: sig
    integer :: ii
    
    select case (Face)
    ! BC on the left face
    case (LEFT)

      select case (dir)
      !BC in X-direction
      case (IAXIS)
        if (xlbc_int == PERIODIC) then
          !do ii  = 1, Gpts
          !  ! lower face
          !  !q(ii, :) = q(iGhi-2*Gpts+ii, :)
          !end do

        else if (xlbc_int == FLOW) then
          do ii  = 1, Gpts
            ! lower face
            q(ilo-ii, :) = q(ilo, :)
          end do

        else if (xlbc_int == REFLECT) then
          sig = 1.000
          if (var == VELX_VAR) sig = -1.000
          do ii = 1, Gpts
            ! lower face
            q(ilo-ii, :) = sig * q(ilo+ii-1, :)
          end do
        end if

      ! BC in Y-direction  
      case (JAXIS)
        if (ylbc_int == PERIODIC) then
          !do ii  = 1, Gpts
          !  ! lower face
          !  !q(:, ii) = q(:, yTpts-2*Gpts+ii)
          !end do

        else if (ylbc_int == FLOW) then
          do ii  = 1, Gpts
            ! lower face
            q(:, jlo-ii) = q(:, jlo)
          end do

        else if (ylbc_int == REFLECT) then
          sig = 1.000
          if (var == VELY_VAR) sig = -1.000
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
    case (RIGHT)
      select case (Dir)
      ! BC in X-direction
      case (IAXIS)
        if (xrbc_int == PERIODIC) then
          !do ii  = 1, Gpts
          !  ! upper face
          !  !q(iGhi-Gpts+ii, :) = q(Gpts + ii, :)
          !end do

        else if (xrbc_int == FLOW) then
          do ii  = 1, Gpts
            ! upper face
            q(ihi+ii, :) = q(ihi, :)
          end do

        else if (xrbc_int == REFLECT) then
          sig = 1.000
          if (var == VELX_VAR) sig = -1.000
          do ii = 1, Gpts
            ! upper face
            q(ihi+ii, :) = sig * q(ihi-ii+1, :)
          end do
        end if

      ! BC in Y-direction
      case (JAXIS)
        if (yrbc_int == PERIODIC) then
          !do ii  = 1, Gpts
          !  ! upper face
          !  !q(:, yTpts-Gpts+ii) = q(:, Gpts + ii)
          !end do

        else if (yrbc_int == FLOW) then
          do ii  = 1, Gpts
            ! upper face
            q(:, jhi+ii) = q(:, jhi)
          end do

        else if (yrbc_int == REFLECT) then
          sig = 1.000
          if (var == VELY_VAR) sig = -1.000
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
    use sim_data, only: mainVar, iGlo, &
                        jGlo, at_xlboundary, &
                        at_ylboundary, at_xrboundary, &
                        at_yrboundary
    implicit none
    real, pointer :: q(:,:)
    integer :: n
    logical :: notappliedxl, notappliedxr, &
               notappliedyl, notappliedyr
    
    ! Periodic BC is already applied in the guardcell_fill routine
    ! applyBC_all should always be called after guardcell_fill routine

    ! defaults for custom BC is false
    notappliedxl=.true.
    notappliedxr=.true.
    notappliedyl=.true.
    notappliedyr=.true.

    if ((at_xlboundary) .and. (notappliedxl)) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          call applyBC(n,q,IAXIS,LEFT)
         nullify(q) 
       end do
    endif

    if ((at_xrboundary) .and. (notappliedxr)) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          call applyBC(n,q,IAXIS,RIGHT)
          nullify(q) 
       end do
    endif

    if ((at_ylboundary) .and. (notappliedyl)) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          call applyBC(n,q,JAXIS,LEFT)
          nullify(q) 
       end do
    endif
    
    if ((at_yrboundary) .and. (notappliedyr)) then 
       do n=NVAR_BEGIN, NVAR_NUMBER
          q(iGlo:,jGlo:) => mainVar(n,:,:)
          call applyBC(n,q,JAXIS,RIGHT)
          nullify(q) 
       end do
    endif

  end subroutine applyBC_all
end module applyBC_module

