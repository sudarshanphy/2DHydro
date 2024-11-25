module recon_module
    implicit none
contains
  function minmod(a,b) result(c)
    implicit none
    real(8) :: a, b
    real(8) :: c, fact

    fact = sign(1.0, a) + sign(1.0, b)
    c = 0.500 * fact * min(abs(a), abs(b)) 

  end function minmod

  subroutine get_facevalue_weno(qm1, q, qp1, fplus12, fminus12, check)
        implicit none
        real(8), intent(in) :: qm1, q, qp1
        real(8), intent(out):: fplus12, fminus12
        logical, intent(in), optional :: check

        real(8) :: epsilon
        real(8), dimension(2) :: beta, alpha, alphabar, &
                                 omega, omegabar, c0, c1, &
                                 c_1, dcoeff, dbarcoeff
        real(8) :: qpluspoly1, qpluspoly2
        real(8) :: qminuspoly1, qminuspoly2
        real(8) :: qplus12, qminus12
        integer :: i
        real(8) :: delp, delm

        dcoeff = (/2.00/3.00, 1.00/3.00/)
        dbarcoeff = (/1.00/3.00, 2.00/3.00/)
        epsilon = 1.0e-6
        beta(1) = (qp1 - q)**2
        beta(2) = (q - qm1)**2

        do i = 1, 2
            alpha(i) = dcoeff(i) / (beta(i) + epsilon)**2
            alphabar(i) = dbarcoeff(i) / (beta(i) + epsilon)**2
        end do

        omega(1) = alpha(1) / sum(alpha)
        omega(2) = alpha(2) / sum(alpha)
        omegabar(1) = alphabar(1) / sum(alphabar)
        omegabar(2) = alphabar(2) / sum(alphabar)

        c0 = (/1.00/2.00, 1.00/2.00/)
        c1 = (/-1.00/2.00, 3.00/2.00/)
        c_1 = (/3.00/2.00, -1.00/2.00/)

        qpluspoly1 = c0(1) * q + c0(2) * qp1
        qpluspoly2 = c1(1) * qm1 + c1(2) * q

        qminuspoly1 = c_1(1) * q + c_1(2) * qp1
        qminuspoly2 = c0(1) * qm1 + c0(2) * q

        qplus12 = omega(1) * qpluspoly1 + omega(2) * qpluspoly2
        qminus12 = omegabar(1) * qminuspoly1 + omegabar(2) * qminuspoly2

        fplus12 = qplus12
        fminus12 = qminus12
        ! monotonicity check
        if ((qplus12 - q) * (q - qminus12) < 0.000) then
           fplus12 = q
           fminus12 = q
        end if
        if (present(check)) then
           delp = qplus12 - q
           delm = q - qminus12
           ! tvd reconstruction
           if (fplus12 < 0.000) then
             fplus12 = q + 0.500 * minmod(delp, delm)
           end if
           if (fminus12 < 0.000) then
             fminus12 = q - 0.500 * minmod(delp, delm)
           end if
        end if   

  end subroutine get_facevalue_weno

  subroutine recon_getcellfaces(dens, momx, momy, ener, dt, &
                                   !method, slimiter, &
                                   xr_plus, xru_plus, xrv_plus, xe_plus, &
                                   xr_minus, xru_minus, xrv_minus, xe_minus, &
                                   yr_plus, yru_plus, yrv_plus, ye_plus, &
                                   yr_minus, yru_minus, yrv_minus, ye_minus)

#include "header.h"        
        use sim_data, only: xTpts, yTpts, dx, dy, gamma, ilo, ihi, jlo, jhi
        implicit none
        real, intent(in) :: dt
        real, dimension(xTpts,yTpts), intent(in) :: dens, momx, momy, ener
        !character(len=*) :: method
        !character(len=*) :: slimiter

        real, dimension(xTpts, yTpts), intent(out) :: xr_plus, xru_plus, xrv_plus, xe_plus
        real, dimension(xTpts, yTpts), intent(out) :: xr_minus, xru_minus, xrv_minus, xe_minus
        real, dimension(xTpts, yTpts), intent(out) :: yr_plus, yru_plus, yrv_plus, ye_plus
        real, dimension(xTpts, yTpts), intent(out) :: yr_minus, yru_minus, yrv_minus, ye_minus

        integer :: i, j

        xr_plus = 0.0
        xru_plus = 0.0
        xrv_plus = 0.0
        xe_plus = 0.0

        xr_minus = 0.0
        xru_minus = 0.0
        xrv_minus = 0.0
        xe_minus = 0.0

        yr_plus = 0.0
        yru_plus = 0.0
        yrv_plus = 0.0
        ye_plus = 0.0

        yr_minus = 0.0
        yru_minus = 0.0
        yrv_minus = 0.0
        ye_minus = 0.0

        do i = ilo-1, ihi+1
            do j = jlo-1, jhi+1

                call get_facevalue_weno(dens(i - 1, j), dens(i, j), dens(i + 1, j), xr_plus(i, j), xr_minus(i, j))
#ifdef DEBUG_WENO
                if ((i == ihi) .and. (j == jlo)) then
                  print *, "hi dens i-1, i, i+1 = ", dens(i-1,j), dens(i,j), dens(i+1,j)
                  print *, "xr_plus, xr_minus = ", xr_plus(i,j), xr_minus(i,j)
                endif
                if ((i == ilo) .and. (j == jlo)) then
                  print *, "lo dens i-1, i, i+1 = ", dens(i-1,j), dens(i,j), dens(i+1,j)
                  print *, "xr_plus, xr_minus = ", xr_plus(i,j), xr_minus(i,j)
                endif
#endif
                call get_facevalue_weno(momx(i - 1, j), momx(i, j), momx(i + 1, j), xru_plus(i, j), xru_minus(i, j))
                call get_facevalue_weno(momy(i - 1, j), momy(i, j), momy(i + 1, j), xrv_plus(i, j), xrv_minus(i, j))
                call get_facevalue_weno(ener(i - 1, j), ener(i, j), ener(i + 1, j), xe_plus(i, j), xe_minus(i, j))

                call get_facevalue_weno(dens(i, j - 1), dens(i, j), dens(i, j + 1), yr_plus(i, j), yr_minus(i, j))
                call get_facevalue_weno(momx(i, j - 1), momx(i, j), momx(i, j + 1), yru_plus(i, j), yru_minus(i, j))
                call get_facevalue_weno(momy(i, j - 1), momy(i, j), momy(i, j + 1), yrv_plus(i, j), yrv_minus(i, j))
                call get_facevalue_weno(ener(i, j - 1), ener(i, j), ener(i, j + 1), ye_plus(i, j), ye_minus(i, j))

            end do
        end do

    end subroutine recon_getcellfaces

end module recon_module
