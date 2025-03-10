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
       use guard_func, only: guardcell_fill
       implicit none
       real, pointer :: solnVar(:,:,:)
       real(8), intent(in) :: dt

       real(8), dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: U, Up1
       real(8), dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: x_plus, x_minus, y_plus, y_minus
                                           
       real(8), dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: xF, yF
       real(8), dimension(NCONSVAR_NUMBER) :: Uleft, Uright, Vleft, Vright
       real(8), dimension(NCONSVAR_NUMBER) :: sgrav, sterm

       integer :: i, j, k, l, m, n
       
       solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)

       sterm(:) = 0.0
       if (usegrav) then
         sterm(3) = grav
         sterm(5) = grav
       end if
       U = 0.0; Up1 = 0.0

       do j = jGlo, jGhi
         do i = iGlo, iGhi
           Up1(DENS_CONS,i,j) = solnVar(DENS_VAR,i,j) 
           Up1(MOMX_CONS:MOMZ_CONS,i,j) = solnVar(DENS_VAR,i,j) * solnVar(VELX_VAR:VELZ_VAR,i,j)
           Up1(ENER_CONS,i,j) = solnVar(ENER_VAR,i,j)

#ifdef MHD
           Up1(BMFX_CONS:BMFZ_CONS,i,j) = solnVar(BMFX_VAR:BMFZ_VAR,i,j)
           Up1(BPSI_CONS,i,j) = solnVar(BPSI_VAR,i,j)
#endif
           U(DENS_CONS,i,j) = solnVar(DENS_VAR,i,j) 
           U(MOMX_CONS:MOMZ_CONS,i,j) = solnVar(DENS_VAR,i,j) * solnVar(VELX_VAR:VELZ_VAR,i,j)
           U(ENER_CONS,i,j) = solnVar(ENER_VAR,i,j)

#ifdef MHD
           U(BMFX_CONS:BMFZ_CONS,i,j) = solnVar(BMFX_VAR:BMFZ_VAR,i,j)
           U(BPSI_CONS,i,j) = solnVar(BPSI_VAR,i,j)
#endif
         end do
       end do
       nullify(solnVar)

       ! rk2 has 2 steps
       do k = 1, 2
         ! use the updated solution for the next step
         solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)

         call recon_getcellfaces(dt, solnVar, &
                                 x_plus, x_minus, y_plus, y_minus)

          do j=jlo, jhi + 1
            do i = ilo, ihi + 1

               do n = 1, NCONSVAR_NUMBER
                 Uleft(n) = x_plus(n,i-1,j)
                 Uright(n) = x_minus(n,i,j)
                 Vleft(n) = y_plus(n,i,j-1)
                 Vright(n) = y_minus(n,i,j)
               end do

              if (to_upper(trim(flux_solver)) == "HLLC") then
#ifndef MHD
                call hllc(Uleft, Uright, "x", xF(:,i,j))
                call hllc(Vleft, Vright, "y", yF(:,i,j))
#else
                print *, "HLLC is only for pure Hydro simulation"
                stop
#endif
              else if (to_upper(trim(flux_solver)) == "HLLE") then
                call hlle(Uleft, Uright, "x", xF(:,i,j))
                call hlle(Vleft, Vright, "y", yF(:,i,j))
              else if (to_upper(trim(flux_solver)) == "HLLD") then
#ifdef MHD
                call hlld(Uleft, Uright, "x", xF(:,i,j))
                call hlld(Vleft, Vright, "y", yF(:,i,j))
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
               sgrav = sterm * (/1.0e0, 0.0, Up1(DENS_CONS,i,j), 0.0, Up1(MOMY_CONS,i,j), &
                                 0.0,0.0,0.0,0.0/)
#else
               sgrav = sterm * (/1.0e0, 0.0, Up1(DENS_CONS,i,j), 0.0, Up1(MOMY_CONS,i,j)/)
#endif
               Up1(:,i,j) = Up1(:,i,j) &
                            + (dt/dx) * (xF(:,i,j) - xF(:,i+1,j)) &
                            + (dt/dy) * (yF(:,i,j) - yF(:,i,j+1)) &
                            + dt * sgrav(:)

               ! get primitive quantities
               solnVar(DENS_VAR,i,j) = Up1(DENS_CONS,i,j)
               solnVar(VELX_VAR:VELZ_VAR,i,j) = Up1(MOMX_CONS:MOMZ_CONS,i,j)/Up1(DENS_CONS,i,j)
               solnVar(ENER_VAR,i,j) = Up1(ENER_CONS,i,j)
#ifdef MHD                        
               solnVar(BMFX_VAR:BMFZ_VAR,i,j) = Up1(BMFX_CONS:BMFZ_CONS,i,j) 
               solnVar(BPSI_VAR,i,j) = Up1(BPSI_CONS,i,j)
#endif
               call eos_getp(solnVar(:,i,j))
            end do
          end do
          nullify(solnVar)
          call guardcell_fill()
          call applyBC_all()

       end do ! rk step loop
       
       solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)

       do j=jGlo, jGhi
         do i = iGlo, iGhi
           Up1(:,i,j) = 0.5e0 * (U(:,i,j) + Up1(:,i,j))

           ! get primitive quantities
           solnVar(DENS_VAR,i,j) = Up1(DENS_CONS,i,j)
           solnVar(VELX_VAR:VELZ_VAR,i,j) = Up1(MOMX_CONS:MOMZ_CONS,i,j)/Up1(DENS_CONS,i,j)
           solnVar(ENER_VAR,i,j) = Up1(ENER_CONS,i,j)
#ifdef MHD                    
           solnVar(BMFX_VAR:BMFZ_VAR,i,j) = Up1(BMFX_CONS:BMFZ_CONS,i,j) 
           solnVar(BPSI_VAR,i,j) = Up1(BPSI_CONS,i,j)
#endif
           call eos_getp(solnVar(:,i,j))

         end do
       end do

       nullify(solnVar) 
    end subroutine RK2_SSP
end module rk2_module
