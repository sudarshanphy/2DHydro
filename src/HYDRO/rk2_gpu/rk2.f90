module rk2_module
#include "param.h"
  implicit none
contains

    subroutine RK2_SSP(dt)
       use flux_module, only: riemann_solve
       use recon_module, only: recon_getcellfaces
       use eos_module, only: eos_getp
       use sim_data, only: grav, usegrav, ilo, ihi, jlo, jhi, &
                           dx, dy, mainVar, iGlo, jGlo, iGhi, &
                           jGhi, smalld, smallp, smalle, gamma

       use applyBC_module, only: applyBC_all 
       use guard_func, only: guardcell_fill
       implicit none
       real, pointer :: solnVar(:,:,:)
       real(8), intent(in) :: dt

       real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: U, Up1
       real(8), dimension(1:NDIM,1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: recon_plus, recon_minus
                                           
       real(8), dimension(1:NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi) :: xF, yF
       real(8), dimension(1:NCONSVAR_NUMBER) :: sgrav, sterm

       integer :: i, j, k, n 
       
       solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)

       sterm(:) = 0.0
       if (usegrav) then
         sterm(3) = grav
         sterm(5) = grav
       end if
       U = 0.0; Up1 = 0.0

       !$OMP TARGET ENTER DATA MAP(ALWAYS, TO: solnVar, dt)
       !$OMP TARGET ENTER DATA MAP(TO: gamma, smalld, smallp, smalle, sterm, dx, dy)
       !$OMP TARGET ENTER DATA MAP(ALLOC: recon_plus, recon_minus, xF, yF, U, Up1, sgrav)
        
       ! Initialize the conserved arrays
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) 
       do j = jGlo, jGhi
         do i = iGlo, iGhi
           Up1(DENS_CONS,i,j) = max(solnVar(DENS_VAR,i,j), smalld) 
           Up1(MOMX_CONS:MOMZ_CONS,i,j) = solnVar(DENS_VAR,i,j) * solnVar(VELX_VAR:VELZ_VAR,i,j)
           Up1(ENER_CONS,i,j) = max(solnVar(ENER_VAR,i,j), smalle)

#ifdef MHD
           Up1(BMFX_CONS:BMFZ_CONS,i,j) = solnVar(BMFX_VAR:BMFZ_VAR,i,j)
           Up1(BPSI_CONS,i,j) = solnVar(BPSI_VAR,i,j)
#endif
           U(DENS_CONS,i,j) = max(solnVar(DENS_VAR,i,j), smalld) 
           U(MOMX_CONS:MOMZ_CONS,i,j) = solnVar(DENS_VAR,i,j) * solnVar(VELX_VAR:VELZ_VAR,i,j)
           U(ENER_CONS,i,j) = max(solnVar(ENER_VAR,i,j), smalle)

#ifdef MHD
           U(BMFX_CONS:BMFZ_CONS,i,j) = solnVar(BMFX_VAR:BMFZ_VAR,i,j)
           U(BPSI_CONS,i,j) = solnVar(BPSI_VAR,i,j)
#endif
         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO 
       nullify(solnVar)

       ! rk2 has 2 steps
       do k = 1, 2
         ! use the updated solution for the next step
         solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)
         if (k == 2) then
            !$OMP TARGET UPDATE TO(solnVar)
         end if

         call recon_getcellfaces(solnVar, &
                                 recon_plus, recon_minus)
         call riemann_solve(recon_plus, recon_minus, xF, yF)

          
          !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
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
               solnVar(DENS_VAR,i,j) = max(Up1(DENS_CONS,i,j), smalld)
               solnVar(VELX_VAR:VELZ_VAR,i,j) = Up1(MOMX_CONS:MOMZ_CONS,i,j)/Up1(DENS_CONS,i,j)
               solnVar(ENER_VAR,i,j) = max(Up1(ENER_CONS,i,j), smalle)
#ifdef MHD                        
               solnVar(BMFX_VAR:BMFZ_VAR,i,j) = Up1(BMFX_CONS:BMFZ_CONS,i,j) 
               solnVar(BPSI_VAR,i,j) = Up1(BPSI_CONS,i,j)
#endif
               ! Apply EOS: in 2D
               solnVar(PRES_VAR,i,j) = (gamma - 1.0) * (solnVar(ENER_VAR,i,j) - 0.50 * solnVar(DENS_VAR,i,j) * &
                                       (solnVar(VELX_VAR,i,j)**2 + solnVar(VELY_VAR,i,j)**2 + &
                                        solnVar(VELZ_VAR,i,j)**2))
#ifdef MHD
               solnVar(PRES_VAR,i,j) = solnVar(PRES_VAR,i,j) - (gamma - 1.0) * 0.50 * &
                                       (solnVar(BMFX_VAR,i,j)**2 + solnVar(BMFY_VAR,i,j)**2 + &
                                        solnVar(BMFZ_VAR,i,j)**2)
#endif
               solnVar(PRES_VAR,i,j) = max(solnVar(PRES_VAR,i,j), smallp)
            end do
          end do
          !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO 
          !$OMP TARGET UPDATE FROM(solnVar)
          nullify(solnVar)
          call guardcell_fill()
          call applyBC_all()

       end do ! rk step loop

       ! get the updated solution and store it as Up1
       solnVar(1:,iGlo:,jGlo:) => mainVar(1:,:,:)

       !$OMP TARGET UPDATE TO(solnVar)
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) 
       do j = jGlo, jGhi
         do i = iGlo, iGhi
           Up1(DENS_CONS,i,j) = max(solnVar(DENS_VAR,i,j), smalld) 
           Up1(MOMX_CONS:MOMZ_CONS,i,j) = solnVar(DENS_VAR,i,j) * solnVar(VELX_VAR:VELZ_VAR,i,j)
           Up1(ENER_CONS,i,j) = max(solnVar(ENER_VAR,i,j), smalle)

#ifdef MHD
           Up1(BMFX_CONS:BMFZ_CONS,i,j) = solnVar(BMFX_VAR:BMFZ_VAR,i,j)
           Up1(BPSI_CONS,i,j) = solnVar(BPSI_VAR,i,j)
#endif
         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO 

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) 
       do j=jGlo, jGhi
         do i = iGlo, iGhi
           Up1(:,i,j) = 0.5e0 * (U(:,i,j) + Up1(:,i,j))

           ! get primitive quantities
           solnVar(DENS_VAR,i,j) = max(Up1(DENS_CONS,i,j), smalld)
           solnVar(VELX_VAR:VELZ_VAR,i,j) = Up1(MOMX_CONS:MOMZ_CONS,i,j)/Up1(DENS_CONS,i,j)
           solnVar(ENER_VAR,i,j) = max(Up1(ENER_CONS,i,j), smalle)
#ifdef MHD                    
           solnVar(BMFX_VAR:BMFZ_VAR,i,j) = Up1(BMFX_CONS:BMFZ_CONS,i,j) 
           solnVar(BPSI_VAR,i,j) = Up1(BPSI_CONS,i,j)
#endif
           ! Apply EOS: in 2D
           solnVar(PRES_VAR,i,j) = (gamma - 1.0) * (solnVar(ENER_VAR,i,j) - 0.50 * solnVar(DENS_VAR,i,j) * &
                                   (solnVar(VELX_VAR,i,j)**2 + solnVar(VELY_VAR,i,j)**2 + &
                                    solnVar(VELZ_VAR,i,j)**2))
#ifdef MHD
           solnVar(PRES_VAR,i,j) = solnVar(PRES_VAR,i,j) - (gamma - 1.0) * 0.50 * &
                                   (solnVar(BMFX_VAR,i,j)**2 + solnVar(BMFY_VAR,i,j)**2 + &
                                    solnVar(BMFZ_VAR,i,j)**2)
#endif
           solnVar(PRES_VAR,i,j) = max(solnVar(PRES_VAR,i,j), smallp)

         end do
       end do
       !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO 
       !$OMP TARGET UPDATE FROM(solnVar)
       !$OMP TARGET EXIT DATA MAP(DELETE: recon_plus, recon_minus, xF, yF)
       
       nullify(solnVar) 
    end subroutine RK2_SSP
end module rk2_module
