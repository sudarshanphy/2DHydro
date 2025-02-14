module rk2_module
#include "header.h"     
#include "param.h"
  implicit none
contains

    subroutine RK2_SSP(dt, solnVar)
#ifndef MHD       
       use riemann_module, only: hllc,hlle
#else
       use riemann_module, only: hlle, hlld
#endif
       use recon_module, only: recon_getcellfaces
       use eos_module, only: eos_getp
       use sim_data, only: grav, usegrav, ilo, ihi, jlo, jhi, &
                           gamma, dx, dy, yTpts, xTpts, smallf, &
                           flux_solver
       use applyBC_module, only: applyBC_all 
       use misc_module, only: to_upper
       implicit none
       real(8), dimension(xTpts, yTpts, NVAR_NUMBER), intent(inout) :: solnVar
       real(8), intent(in) :: dt

       real(8), dimension(xTpts, yTpts, NCONSVAR_NUMBER) :: U, Up1
       real(8), dimension(xTpts, yTpts, NVAR_NUMBER) :: x_plus, x_minus, y_plus, y_minus
   
       real(8), dimension(xTpts, yTpts, NCONSVAR_NUMBER) :: xF, yF
       real(8), dimension(NVAR_NUMBER-1) :: Uleft, Uright, Vleft, Vright
       real(8), dimension(NCONSVAR_NUMBER) :: sgrav, sterm

       integer :: i, j, k, l, m, n

#ifdef DEBUG_PAR
       print *, "usegrav, grav = ", usegrav, grav
       print *, "ilo, ihi, jlo, jhi = ", ilo, ihi, jlo, jhi
       print *, "xTpts, yTpts, dx, dy = ", xTpts, yTpts, dx, dy
       print *, "gamma = ", gamma
#endif

       sterm(:) = 0.0
       if (usegrav) then
         sterm(3) = grav
         sterm(5) = grav
       end if
       U = 0.0; Up1 = 0.0

       do j = 1, yTpts
         do i = 1, xTpts
           U(i,j,DENS_CONS) = solnVar(i,j,DENS_VAR) 
           U(i,j,MOMX_CONS) = solnVar(i,j,DENS_VAR) * solnVar(i,j,VELX_VAR)
           U(i,j,MOMY_CONS) = solnVar(i,j,DENS_VAR) * solnVar(i,j,VELY_VAR)
           U(i,j,MOMZ_CONS) = solnVar(i,j,DENS_VAR) * solnVar(i,j,VELZ_VAR)
           U(i,j,ENER_CONS) = solnVar(i,j,ENER_VAR)

#ifdef MHD
           U(i,j,BMFX_CONS) = solnVar(i,j,BMFX_VAR)
           U(i,j,BMFY_CONS) = solnVar(i,j,BMFY_VAR)
           U(i,j,BMFZ_CONS) = solnVar(i,j,BMFZ_VAR)
           U(i,j,BPSI_CONS) = solnVar(i,j,BPSI_VAR)
#endif
         end do
       end do
      
       ! rk2 has 2 steps
       do k = 1, 2

         call recon_getcellfaces(dt, solnVar, &
                                 x_plus, x_minus, y_plus, y_minus)

          ! initialize them to 0
          xF(:,:,:) = 0.0
          yF(:,:,:) = 0.0

          do j=jlo, jhi + 1
            do i = ilo, ihi + 1

               do n = 1, NVAR_NUMBER - 1
                 Uleft(n) = x_plus(i-1,j,n)
                 Uright(n) = x_minus(i,j,n)
                 Vleft(n) = y_plus(i,j-1,n)
                 Vright(n) = y_minus(i,j,n)

               end do

              if (to_upper(trim(flux_solver)) == "HLLC") then
#ifndef MHD
                call hllc(Uleft, Uright, "x", xF(i,j,:))
                call hllc(Vleft, Vright, "y", yF(i,j,:))
#else
                print *, "HLLC is only for pure Hydro simulation"
                stop
#endif
              else if (to_upper(trim(flux_solver)) == "HLLE") then
                call hlle(Uleft, Uright, "x", xF(i,j,:))
                call hlle(Vleft, Vright, "y", yF(i,j,:))
              else if (to_upper(trim(flux_solver)) == "HLLD") then
#ifdef MHD
                call hlld(Uleft, Uright, "x", xF(i,j,:))
                call hlld(Vleft, Vright, "y", yF(i,j,:))
#else
                print *, "HLLD is only for MHD simulation"
                stop
#endif
              else
                print *, "Available solvers:"
                print *, "HLLE/HLLC -> Hydro"
                print *, "HLLE/HLLD -> MHD"
                stop
              end if
             
               
            end do
          end do

          do j=jlo, jhi
            do i = ilo, ihi
               do n = 1, NCONSVAR_NUMBER
                  !sgrav = sterm * (/1.0e0, dens(i,j), dens(i,j), 0.0, momy(i,j),0.0,0.0,0.0,0.0/)
                  Up1(i,j,n) = U(i,j,n) &
                               + (dt/dx) * (xF(i,j,n) - xF(i+1,j,n)) &
                               + (dt/dy) * (yF(i,j,n) - yF(i,j+1,n))
               end do
!
               solnVar(i,j,DENS_VAR) = Up1(i,j,DENS_CONS)
               solnVar(i,j,VELX_VAR) = Up1(i,j,MOMX_CONS)/Up1(i,j,DENS_CONS)
               solnVar(i,j,VELY_VAR) = Up1(i,j,MOMY_CONS)/Up1(i,j,DENS_CONS)
               solnVar(i,j,ENER_VAR) = Up1(i,j,ENER_CONS)
#ifdef MHD                        
               solnVar(i,j,BMFX_VAR) = Up1(i,j,BMFX_CONS) 
               solnVar(i,j,BMFY_VAR) = Up1(i,j,BMFY_CONS)
               solnVar(i,j,BPSI_VAR) = Up1(i,j,BMFZ_CONS)
#endif
               solnVar(i,j,PRES_VAR) = eos_getp(solnVar(i,j,:)) 
            end do
          end do
          call applyBC_all(solnVar)
       end do

       do j=jlo, jhi
         do i = ilo, ihi
           do n = 1, NCONSVAR_NUMBER
              Up1(i,j,n) = 0.5e0 * (U(i,j,n) + Up1(i,j,n))
           end do
           ! get primitive quantities
           solnVar(i,j,DENS_VAR) = Up1(i,j,DENS_CONS)
           solnVar(i,j,VELX_VAR) = Up1(i,j,MOMX_CONS)/Up1(i,j,DENS_CONS)
           solnVar(i,j,VELY_VAR) = Up1(i,j,MOMY_CONS)/Up1(i,j,DENS_CONS)
           solnVar(i,j,ENER_VAR) = Up1(i,j,ENER_CONS)
#ifdef MHD                    
           solnVar(i,j,BMFX_VAR) = Up1(i,j,BMFX_CONS) 
           solnVar(i,j,BMFY_VAR) = Up1(i,j,BMFY_CONS)
           solnVar(i,j,BPSI_VAR) = Up1(i,j,BMFZ_CONS)
#endif
           solnVar(i,j,PRES_VAR) = eos_getp(solnVar(i,j,:)) 
         end do
       end do
       
    end subroutine RK2_SSP
end module rk2_module
