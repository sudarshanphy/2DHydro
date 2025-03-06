module rk2_module
#include "param.h"
  implicit none
contains

    subroutine RK2_SSP(dt)
#ifndef MHD       
       use flux_module, only: hllc,hlle
#else
       use flux_module, only: hlle, hlld
#endif
       use recon_module, only: recon_getcellfaces
       use eos_module, only: eos_getp
       use sim_data, only: grav, usegrav, ilo, ihi, jlo, jhi, &
                           gamma, dx, dy, lyTpts, lxTpts, smallf, &
                           flux_solver, mainVar, iGlo, jGlo, iGhi, &
                           jGhi

       use applyBC_module, only: applyBC_all 
       use misc_module, only: to_upper
       implicit none
       real, pointer :: solnVar(:,:,:)
       real(8), intent(in) :: dt

       real(8), dimension(iGlo:iGhi, jGlo:jGhi, NCONSVAR_NUMBER) :: U, Up1
       real(8), dimension(iGlo:iGhi,jGlo:jGhi, NCONSVAR_NUMBER) :: x_plus, x_minus, y_plus, y_minus
   
       real(8), dimension(iGlo:iGhi, jGlo:jGhi, NCONSVAR_NUMBER) :: xF, yF
       real(8), dimension(NCONSVAR_NUMBER) :: Uleft, Uright, Vleft, Vright
       real(8), dimension(NCONSVAR_NUMBER) :: sgrav, sterm

       integer :: i, j, k, l, m, n
       
       solnVar(iGlo:,jGlo:,1:) => mainVar(:,:,1:)

       sterm(:) = 0.0
       if (usegrav) then
         sterm(3) = grav
         sterm(5) = grav
       end if
       U = 0.0; Up1 = 0.0

       do j = jGlo, jGhi
         do i = iGlo, iGhi
           Up1(i,j,DENS_CONS) = solnVar(i,j,DENS_VAR) 
           Up1(i,j,MOMX_CONS:MOMZ_CONS) = solnVar(i,j,DENS_VAR) * solnVar(i,j,VELX_VAR:VELZ_VAR)
           Up1(i,j,ENER_CONS) = solnVar(i,j,ENER_VAR)

#ifdef MHD
           Up1(i,j,BMFX_CONS:BMFZ_CONS) = solnVar(i,j,BMFX_VAR:BMFZ_VAR)
           Up1(i,j,BPSI_CONS) = solnVar(i,j,BPSI_VAR)
#endif
           U(i,j,DENS_CONS) = solnVar(i,j,DENS_VAR) 
           U(i,j,MOMX_CONS:MOMZ_CONS) = solnVar(i,j,DENS_VAR) * solnVar(i,j,VELX_VAR:VELZ_VAR)
           U(i,j,ENER_CONS) = solnVar(i,j,ENER_VAR)

#ifdef MHD
           U(i,j,BMFX_CONS:BMFZ_CONS) = solnVar(i,j,BMFX_VAR:BMFZ_VAR)
           U(i,j,BPSI_CONS) = solnVar(i,j,BPSI_VAR)
#endif
         end do
       end do
       nullify(solnVar)

       ! rk2 has 2 steps
       do k = 1, 2
         ! use the updated solution for the next step
         solnVar(iGlo:,jGlo:,1:) => mainVar(:,:,1:)

         call recon_getcellfaces(dt, solnVar, &
                                 x_plus, x_minus, y_plus, y_minus)

          do j=jlo, jhi + 1
            do i = ilo, ihi + 1

               do n = 1, NCONSVAR_NUMBER
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
#ifdef MHD
               sgrav = sterm * (/1.0e0, 0.0, Up1(i,j,DENS_CONS), 0.0, Up1(i,j,MOMY_CONS), &
                                 0.0,0.0,0.0,0.0/)
#else
               sgrav = sterm * (/1.0e0, 0.0, Up1(i,j,DENS_CONS), 0.0, Up1(i,j,MOMY_CONS)/)
#endif
               Up1(i,j,:) = Up1(i,j,:) &
                            + (dt/dx) * (xF(i,j,:) - xF(i+1,j,:)) &
                            + (dt/dy) * (yF(i,j,:) - yF(i,j+1,:)) &
                            + dt * sgrav(:)

               ! get primitive quantities
               solnVar(i,j,DENS_VAR) = Up1(i,j,DENS_CONS)
               solnVar(i,j,VELX_VAR:VELZ_VAR) = Up1(i,j,MOMX_CONS:MOMZ_CONS)/Up1(i,j,DENS_CONS)
               solnVar(i,j,ENER_VAR) = Up1(i,j,ENER_CONS)
#ifdef MHD                        
               solnVar(i,j,BMFX_VAR:BMFZ_VAR) = Up1(i,j,BMFX_CONS:BMFZ_CONS) 
               solnVar(i,j,BPSI_VAR) = Up1(i,j,BPSI_CONS)
#endif
               call eos_getp(solnVar(i,j,:))
            end do
          end do
          nullify(solnVar)
          call guardcell_fill()
          call applyBC_all()

       end do ! rk step loop
       
       solnVar(iGlo:,jGlo:,1:) => mainVar(:,:,1:)

       do j=jGlo, jGhi
         do i = iGlo, iGhi
           Up1(i,j,:) = 0.5e0 * (U(i,j,:) + Up1(i,j,:))

           ! get primitive quantities
           solnVar(i,j,DENS_VAR) = Up1(i,j,DENS_CONS)
           solnVar(i,j,VELX_VAR:VELZ_VAR) = Up1(i,j,MOMX_CONS:MOMZ_CONS)/Up1(i,j,DENS_CONS)
           solnVar(i,j,ENER_VAR) = Up1(i,j,ENER_CONS)
#ifdef MHD                    
           solnVar(i,j,BMFX_VAR:BMFZ_VAR) = Up1(i,j,BMFX_CONS:BMFZ_CONS) 
           solnVar(i,j,BPSI_VAR) = Up1(i,j,BPSI_CONS)
#endif
           call eos_getp(solnVar(i,j,:))

         end do
       end do

       nullify(solnVar) 
    end subroutine RK2_SSP
end module rk2_module
