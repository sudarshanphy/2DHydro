module rk2_module
  implicit none
contains

    subroutine RK2_SSP(dens, velx, vely, pres, ener, dt)
#include "header.h"      
       
       use hllc_module, only: hllc
       use recon_module, only: recon_getcellfaces
       use eos_module, only: eos_getp
       use sim_data, only: grav, usegrav, ilo, ihi, jlo, jhi, &
                           gamma, dx, dy, yTpts, xTpts, smallf
       use applyBC_module, only: applyBC_all 
       implicit none
       real(8), dimension(xTpts, yTpts), intent(inout) :: dens, velx, vely, pres, ener
       real(8), intent(in) :: dt 
       real(8), dimension(xTpts, yTpts) :: momx, momy, dens_n, &
                                           momx_n, momy_n, ener_n
       real(8), dimension(xTpts, yTpts) :: xr_plus, xru_plus, xrv_plus, xe_plus
       real(8), dimension(xTpts, yTpts) :: xr_minus, xru_minus, xrv_minus, xe_minus
       real(8), dimension(xTpts, yTpts) :: yr_plus, yru_plus, yrv_plus, ye_plus
       real(8), dimension(xTpts, yTpts) :: yr_minus, yru_minus, yrv_minus, ye_minus

       real(8), dimension(xTpts, yTpts) :: xrF, xruF, xrvF, xrwF, xeF
       real(8), dimension(xTpts, yTpts) :: yrF, yruF, yrvF, yrwF, yeF

       real(8), dimension(5) :: Uleft, Uright, Vleft, Vright, xF, yF
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
             momx(i,j) = velx(i,j) * dens(i,j)
             momy(i,j) = vely(i,j) * dens(i,j)
             dens_n(i,j) = dens(i,j)
             momx_n(i,j) = momx(i,j)
             momy_n(i,j) = momy(i,j)
             ener_n(i,j) = ener(i,j)
         end do
       end do
      
       ! rk2 has 2 steps
       do k = 1, 2

         call recon_getcellfaces(dens, momx, momy, ener, dt, &
                                 !method, slimiter, &
                                 xr_plus, xru_plus, xrv_plus, xe_plus, &
                                 xr_minus, xru_minus, xrv_minus, xe_minus, &
                                 yr_plus, yru_plus, yrv_plus, ye_plus, &
                                 yr_minus, yru_minus, yrv_minus, ye_minus)

          do j=jlo, jhi + 1
            do i = ilo, ihi + 1
               ! initialize them to 0
               xF(:) = 0.0
               yF(:) = 0.0

               Uleft = (/xr_plus(i-1,j), xru_plus(i-1,j), &
                 xrv_plus(i-1,j), 0.0e0, xe_plus(i-1,j)/)
               Uright = (/xr_minus(i,j), xru_minus(i,j), &
                 xrv_minus(i,j), 0.0e0, xe_minus(i,j)/)
               Vleft = (/yr_plus(i,j-1), yru_plus(i,j-1), &
                 yrv_plus(i,j-1), 0.0e0, ye_plus(i,j-1)/)
               Vright = (/yr_minus(i,j), yru_minus(i,j), &
                 yrv_minus(i,j), 0.0e0, ye_minus(i,j)/)

              call hllc(Uleft, Uright, "x", xF)
              call hllc(Vleft, Vright, "y", yF)
              
              xrF(i,j) = xF(1); xruF(i,j) = xF(2); xrvF(i,j) = xF(3); xeF(i,j) = xF(5)
              yrF(i,j) = yF(1); yruF(i,j) = yF(2); yrvF(i,j) = yF(3); yeF(i,j) = yF(5)
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
                         
               ! floor values for minimum is smallf          
               !dens(i,j) = max(dens(i,j), smallf)
               !momx(i,j) = max(abs(momx(i,j)), smallf)
               !momy(i,j) = max(abs(momy(i,j)), smallf)
               !ener(i,j) = max(ener(i,j), smallf)

            end do
          end do

          call applyBC_all(dens, momx, momy, pres, ener)

       end do

       do j=jlo, jhi
         do i = ilo, ihi
           dens(i,j) = 0.5e0 * (dens_n(i,j) + dens(i,j))
           momx(i,j) = 0.5e0 * (momx_n(i,j) + momx(i,j))
           momy(i,j) = 0.5e0 * (momy_n(i,j) + momy(i,j))
           ener(i,j) = 0.5e0 * (ener_n(i,j) + ener(i,j))
         end do
       end do
       
       do j = jlo, jhi
         do i = ilo, ihi
             velx(i,j) = momx(i,j) / dens(i,j)
             vely(i,j) = momy(i,j) / dens(i,j)
             pres(i,j) = eos_getp((/dens(i,j), velx(i,j), vely(i,j), 0.0, ener(i,j)/))

             !velx(i,j) = max(abs(velx(i,j)), smallf)
             !vely(i,j) = max(abs(vely(i,j)), smallf)
             !pres(i,j) = max(pres(i,j), smallf)
         end do
       end do

    end subroutine RK2_SSP
end module rk2_module
