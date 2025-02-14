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
       real(8), dimension(xTpts, yTpts) :: dens, pres, velx, vely, ener 
       real(8), dimension(xTpts, yTpts) :: momx, momy, dens_n, &
                                           momx_n, momy_n, ener_n
       real(8), dimension(xTpts, yTpts) :: xr_plus, xu_plus, xv_plus, xp_plus
       real(8), dimension(xTpts, yTpts) :: xr_minus, xu_minus, xv_minus, xp_minus
       real(8), dimension(xTpts, yTpts) :: yr_plus, yu_plus, yv_plus, yp_plus
       real(8), dimension(xTpts, yTpts) :: yr_minus, yu_minus, yv_minus, yp_minus
#ifdef MHD
       real(8), dimension(xTpts, yTpts) :: bmfx, bmfy, bpsi
       real(8), dimension(xTpts, yTpts) :: bmfx_n, bmfy_n, bpsi_n
       real(8), dimension(xTpts, yTpts) :: xbx_plus, xby_plus, xbx_minus, xby_minus, &
                                           xbp_plus, xbp_minus
       real(8), dimension(xTpts, yTpts) :: ybx_plus, yby_plus, ybx_minus, yby_minus, &
                                           ybp_plus, ybp_minus
#endif

       real(8), dimension(xTpts, yTpts) :: xrF, xruF, xrvF, xrwF, xeF
       real(8), dimension(xTpts, yTpts) :: yrF, yruF, yrvF, yrwF, yeF
#ifdef MHD
       real(8), dimension(xTpts, yTpts) :: xbxF, xbyF, xbpF
       real(8), dimension(xTpts, yTpts) :: ybxF, ybyF, ybpF
#endif

       real(8), dimension(NVAR_NUMBER-1) :: Uleft, Uright, Vleft, Vright, xF, yF
       real(8), dimension(4) :: sgrav, sterm

       integer :: i, j, k, l, m, n

#ifdef DEBUG_PAR
       print *, "usegrav, grav = ", usegrav, grav
       print *, "ilo, ihi, jlo, jhi = ", ilo, ihi, jlo, jhi
       print *, "xTpts, yTpts, dx, dy = ", xTpts, yTpts, dx, dy
       print *, "gamma = ", gamma
#endif

       if (usegrav) then
         sterm(:) = (/0.0e0, 0.0e0, grav, grav/)
       else
         sterm(:) = 0.0e0
       end if
       do j = 1, yTpts
         do i = 1, xTpts
         dens(i,j) = solnVar(i,j,DENS_VAR)
         velx(i,j) = solnVar(i,j,VELX_VAR)
         vely(i,j) = solnVar(i,j,VELY_VAR)
         pres(i,j) = solnVar(i,j,PRES_VAR)
         ener(i,j) = solnVar(i,j,ENER_VAR)
   
#ifdef MHD
       bmfx = solnVar(:,:,BMFX_VAR)
       bmfy = solnVar(:,:,BMFY_VAR)
       bpsi = solnVar(:,:,BPSI_VAR)
#endif
             momx(i,j) = velx(i,j) * dens(i,j)
             momy(i,j) = vely(i,j) * dens(i,j)
             dens_n(i,j) = dens(i,j)
             momx_n(i,j) = momx(i,j)
             momy_n(i,j) = momy(i,j)
             ener_n(i,j) = ener(i,j)
#ifdef MHD
             bmfx_n(i,j) = bmfx(i,j)
             bmfy_n(i,j) = bmfy(i,j)
             bpsi_n(i,j) = bpsi(i,j)
#endif
         end do
       end do
      
       ! rk2 has 2 steps
       do k = 1, 2

       dens(:,:) = solnVar(:,:,DENS_VAR)
       velx(:,:) = solnVar(:,:,VELX_VAR)
       vely(:,:) = solnVar(:,:,VELY_VAR)
       pres(:,:) = solnVar(:,:,PRES_VAR)
       ener(:,:) = solnVar(:,:,ENER_VAR)

#ifdef MHD
       bmfx = solnVar(:,:,BMFX_VAR)
       bmfy = solnVar(:,:,BMFY_VAR)
       bpsi = solnVar(:,:,BPSI_VAR)
#endif

         call recon_getcellfaces(dt, dens, velx, vely, pres, &
#ifdef MHD 
                                           bmfx, bmfy, bpsi, &
#endif
                                 xr_plus, xu_plus, xv_plus, xp_plus, &
#ifdef MHD
                                 xbx_plus, xby_plus, xbp_plus, &
#endif
                                 xr_minus, xu_minus, xv_minus, xp_minus, &
#ifdef MHD
                                 xbx_minus, xby_minus, xbp_minus, &
#endif
                                 yr_plus, yu_plus, yv_plus, yp_plus, &
#ifdef MHD
                                 ybx_plus, yby_plus, ybp_plus, &
#endif
                                 yr_minus, yu_minus, yv_minus, yp_minus &
#ifdef MHD
                                 , ybx_minus, yby_minus, ybp_minus      &
#endif
                                                                        &)

          do j=jlo, jhi + 1
            do i = ilo, ihi + 1
               ! initialize them to 0
               xF(:) = 0.0
               yF(:) = 0.0

               Uleft = (/xr_plus(i-1,j), xu_plus(i-1,j), &
                 xv_plus(i-1,j), 0.0e0, xp_plus(i-1,j)   &
#ifdef MHD
                 , xbx_plus(i-1,j), xby_plus(i-1,j), 0.0, xbp_plus(i-1,j) &
#endif
                                                        &/)
               Uright = (/xr_minus(i,j), xu_minus(i,j), &
                 xv_minus(i,j), 0.0e0, xp_minus(i,j)    &
#ifdef MHD
                 , xbx_minus(i,j), xby_minus(i,j), 0.0, xbp_minus(i,j) &
#endif
                                                         &/)
               Vleft = (/yr_plus(i,j-1), yu_plus(i,j-1), &
                 yv_plus(i,j-1), 0.0e0, yp_plus(i,j-1)   &
#ifdef MHD
                 , ybx_plus(i,j-1), yby_plus(i,j-1), 0.0, ybp_plus(i,j-1) &
#endif
                                                         &/)
               Vright = (/yr_minus(i,j), yu_minus(i,j), &
                 yv_minus(i,j), 0.0e0, yp_minus(i,j)    &
#ifdef MHD
                 , ybx_minus(i,j), yby_minus(i,j), 0.0, ybp_minus(i,j) &
#endif
                                                        &/)

              if (to_upper(trim(flux_solver)) == "HLLC") then
#ifndef MHD
                call hllc(Uleft, Uright, "x", xF)
                call hllc(Vleft, Vright, "y", yF)
#else
                print *, "HLLC is only for pure Hydro simulation"
                stop
#endif
              else if (to_upper(trim(flux_solver)) == "HLLE") then
                call hlle(Uleft, Uright, "x", xF)
                call hlle(Vleft, Vright, "y", yF)
              else if (to_upper(trim(flux_solver)) == "HLLD") then
#ifdef MHD
                call hlld(Uleft, Uright, "x", xF)
                call hlld(Vleft, Vright, "y", yF)
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
              
              xrF(i,j) = xF(1); xruF(i,j) = xF(2); xrvF(i,j) = xF(3); xeF(i,j) = xF(5)
              yrF(i,j) = yF(1); yruF(i,j) = yF(2); yrvF(i,j) = yF(3); yeF(i,j) = yF(5)
#ifdef MHD
              xbxF(i,j) = xF(6); xbyF(i,j) = xF(7); xbpF(i,j) = xF(9)
              ybxF(i,j) = yF(6); ybyF(i,j) = yF(7); ybpF(i,j) = yF(9)
#endif
            end do
          end do

          do j=jlo, jhi
            do i = ilo, ihi
               sgrav = sterm * (/1.0e0, dens(i,j), dens(i,j), momy(i,j)/)
               dens(i,j) = dens(i,j) &
                           + (dt / dx) * (xrF(i,j) - xrF(i+1, j)) &
                           + (dt / dy) * (yrF(i,j) - yrF(i, j+1)) &
                           + dt * sgrav(1)
               momx(i,j) = momx(i,j) &
                           + (dt / dx) * (xruF(i,j) - xruF(i+1, j)) &
                           + (dt / dy) * (yruF(i,j) - yruF(i, j+1)) &
                           + dt * sgrav(2)
               momy(i,j) = momy(i,j) &
                           + (dt / dx) * (xrvF(i,j) - xrvF(i+1, j)) &
                           + (dt / dy) * (yrvF(i,j) - yrvF(i, j+1)) &
                           + dt * sgrav(3)
               ener(i,j) = ener(i,j) &
                           + (dt / dx) * (xeF(i,j) - xeF(i+1, j)) &
                           + (dt / dy) * (yeF(i,j) - yeF(i, j+1)) &
                           + dt * sgrav(4)

#ifdef MHD
               bmfx(i,j) = bmfx(i,j) & 
                           + (dt / dx) * (xbxF(i,j) - xbxF(i+1, j)) &
                           + (dt / dy) * (ybxF(i,j) - ybxF(i, j+1))
               bmfy(i,j) = bmfy(i,j) &
                           + (dt / dx) * (xbyF(i,j) - xbyF(i+1, j)) &
                           + (dt / dy) * (ybyF(i,j) - ybyF(i, j+1))
               bpsi(i,j) = bpsi(i,j) &
                           + (dt / dx) * (xbpF(i,j) - xbpF(i+1, j)) &
                           + (dt / dy) * (ybpF(i,j) - ybpF(i, j+1))
#endif

               velx(i,j) = momx(i,j)/dens(i,j)
               vely(i,j) = momy(i,j)/dens(i,j)
#ifdef MHD
               pres(i,j) = eos_getp((/dens(i,j), velx(i,j), vely(i,j), 0.0, ener(i,j), &
                                      bmfx(i,j), bmfy(i,j), 0.0/))      
#else
               pres(i,j) = eos_getp((/dens(i,j), velx(i,j), vely(i,j), 0.0, ener(i,j)/))      
#endif
          solnVar(i,j,DENS_VAR) = dens(i,j)
          solnVar(i,j,VELX_VAR) = velx(i,j)
          solnVar(i,j,VELY_VAR) = vely(i,j)
          solnVar(i,j,PRES_VAR) = pres(i,j)
          solnVar(i,j,ENER_VAR) = ener(i,j)
                                      
#ifdef MHD                        
          solnVar(i,j,BMFX_VAR) = bmfx(i,j)
          solnVar(i,j,BMFY_VAR) = bmfy(i,j)
          solnVar(i,j,BPSI_VAR) = bpsi(i,j)
#endif
            end do
          end do
          call applyBC_all(solnVar)
       end do

       do j=jlo, jhi
         do i = ilo, ihi
           dens(i,j) = 0.5e0 * (dens_n(i,j) + dens(i,j))
           momx(i,j) = 0.5e0 * (momx_n(i,j) + momx(i,j))
           momy(i,j) = 0.5e0 * (momy_n(i,j) + momy(i,j))
           ener(i,j) = 0.5e0 * (ener_n(i,j) + ener(i,j))

#ifdef MHD
          bmfx(i,j) = 0.5e0 * (bmfx_n(i,j) + bmfx(i,j))
          bmfy(i,j) = 0.5e0 * (bmfy_n(i,j) + bmfy(i,j))
          bpsi(i,j) = 0.5e0 * (bpsi_n(i,j) + bpsi(i,j))
#endif
           ! get primitive quantities
           velx(i,j) = momx(i,j) / dens(i,j)
           vely(i,j) = momy(i,j) / dens(i,j)
#ifdef MHD
           pres(i,j) = eos_getp((/dens(i,j), velx(i,j), vely(i,j), 0.0, ener(i,j), &
                                  bmfx(i,j), bmfy(i,j), 0.0/))      
#else
           pres(i,j) = eos_getp((/dens(i,j), velx(i,j), vely(i,j), 0.0, ener(i,j)/))      
#endif

           solnVar(i,j,DENS_VAR) = dens(i,j)
           solnVar(i,j,VELX_VAR) = velx(i,j)
           solnVar(i,j,VELY_VAR) = vely(i,j)
           solnVar(i,j,PRES_VAR) = pres(i,j)
           solnVar(i,j,ENER_VAR) = ener(i,j)
#ifdef MHD                     
           solnVar(i,j,BMFX_VAR) = bmfx(i,j)
           solnVar(i,j,BMFY_VAR) = bmfy(i,j)
           solnVar(i,j,BPSI_VAR) = bpsi(i,j)
#endif
         end do
       end do
       
    end subroutine RK2_SSP
end module rk2_module
