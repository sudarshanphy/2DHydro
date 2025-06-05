module recon_module
#include "param.h"
    implicit none
contains
  function minmod(a,b) result(c)
    implicit none
    real :: a, b
    real :: c, fact

    fact = sign(1.0, a) + sign(1.0, b)
    c = 0.500 * fact * min(abs(a), abs(b)) 

  end function minmod

  subroutine get_facevalue_weno(solnVar, flxplus, flxminus, dir, check)
        use sim_data, only: ilo, ihi, jlo, jhi, iGlo, iGhi, jGlo, jGhi 
        implicit none
        real, dimension(1:, iGlo:, jGlo:), intent(in) :: solnVar
        real, dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi), intent(inout) :: flxminus, flxplus
        integer, intent(in) :: dir
        logical, intent(in), optional :: check
        real :: fplus12, fminus12

        real, parameter :: epsilon = 1.0e-6
        real, dimension(2) :: beta, alpha, alphabar, &
                                 omega, omegabar, & 
                                 qpluspoly, qminuspoly
        real :: qplus12, qminus12
        integer :: i, j, n, im, ip, jm, jp
        real :: delp, delm
        real :: qm1, q, qp1
        ! coefficients for the polynomial
        real, dimension(2) :: c0 = (/1.00/2.00, 1.00/2.00/)
        real, dimension(2) :: c1 = (/-1.00/2.00, 3.00/2.00/)
        real, dimension(2) :: c_1 = (/3.00/2.00, -1.00/2.00/)
        ! linear weights 
        real, dimension(2) :: dcoeff = (/2.00/3.00, 1.00/3.00/)
        real, dimension(2) :: dbarcoeff = (/1.00/3.00, 2.00/3.00/)

        select case (dir)
           case (1)
              im = 1; ip = 1; jm = 0; jp = 0
           case(2)
              im = 0; ip = 0; jm = 1; jp = 1
        end select

        do n = 1, NCONSVAR_NUMBER
           do j = jlo-1, jhi+1
               do i = ilo-1, ihi+1

                  qm1 = solnVar(n,i-im, j-jm)
                  q = solnVar(n,i,j)
                  qp1 = solnVar(n,i+ip, j+jp)

                  beta(1) = (qp1 - q)**2
                  beta(2) = (q - qm1)**2

                  alpha(1:2) = dcoeff(1:2) / (beta(1:2) + epsilon)**2
                  alphabar(1:2) = dbarcoeff(1:2) / (beta(1:2) + epsilon)**2

                  omega(1:2) = alpha(1:2) / sum(alpha)
                  omegabar(1:2) = alphabar(1:2) / sum(alphabar)

                  qpluspoly(1) = dot_product(c0(1:2), (/q, qp1/))
                  qpluspoly(2) = dot_product(c1(1:2), (/qm1, q/))

                  qminuspoly(1) = dot_product(c_1(1:2), (/q, qp1/))
                  qminuspoly(2) = dot_product(c0(1:2), (/qm1, q/))

                  qplus12 = dot_product(omega(1:2), qpluspoly(1:2))
                  qminus12 = dot_product(omegabar(1:2), qminuspoly(1:2))

                  fplus12 = qplus12
                  fminus12 = qminus12
                  ! monotonicity check
                  if ((qplus12 - q) * (q - qminus12) < 0.000) then
                     fplus12 = q
                     fminus12 = q
                  end if
                  if (present(check)) then
                     delp = qp1 - q
                     delm = q - qm1
                     ! tvd reconstruction
                     if (fplus12 < 0.000) then
                       fplus12 = q + 0.500 * minmod(delp, delm)
                     end if
                     if (fminus12 < 0.000) then
                       fminus12 = q - 0.500 * minmod(delp, delm)
                     end if
                  end if   
                  flxminus(n,i,j) = fminus12
                  flxplus(n,i,j) = fplus12   
               end do
           end do
        end do

  end subroutine get_facevalue_weno

  subroutine get_facevalue_weno5z(solnVar, flxplus, flxminus, dir, check)
        use sim_data, only: ilo, ihi, jlo, jhi, iGlo, iGhi, jGlo, jGhi 
        implicit none
        real, dimension(1:, iGlo:, jGlo:), intent(in) :: solnVar
        real, dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi), intent(inout) :: flxminus, flxplus
        integer, intent(in) :: dir
        logical, intent(in), optional :: check

        real, dimension(1:5) :: qin
        real :: fplus12, fminus12
        real, parameter :: epsilon = 1.0e-36
        real, dimension(3) :: beta, alpha, alphabar, &
                              omega, omegabar, qpluspoly, &
                              qminuspoly
        real :: qplus12, qminus12
        integer :: i, j, n, imm, im, ip, ipp, &
                            jmm, jm, jp, jpp  
        real :: delp, delm, absbeta_diff
        real :: qm2, qm1, q, qp1, qp2

        ! adapted from flashx spark
        real :: n13o12 = 13.0/12.0
        !! Set WENO5 coefficients once and for all
        !u_{1,i+1/2}= 2/6*u_{i-2} -7/6*u_{i-1} +11/6*u_{i}
        real, dimension(3) :: coeff1p1 = (/ 2./6., -7./6., 11./6./)
        !u_{2,i+1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
        real, dimension(3) :: coeff1p2 = (/-1./6.,  5./6.,  2./6./)
        !u_{3,i+1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
        real, dimension(3) :: coeff1p3 = (/ 2./6.,  5./6., -1./6./)
        ! linear weights for i + 1/2
        real, dimension(3) :: dcoeff = (/0.1, 0.6, 0.3/)
        !u_{1,i-1/2}=-1/6*u_{i-2} +5/6*u_{i-1} + 2/6*u_{i}
        real, dimension(3) :: coeff1m1 = (/-1./6.,  5./6.,  2./6./)
        !u_{2,i-1/2}= 2/6*u_{i-2} +5/6*u_{i-1} - 1/6*u_{i}
        real, dimension(3) :: coeff1m2 = (/ 2./6.,  5./6., -1./6./)
        !u_{3,i-1/2}=11/6*u_{i-2} -7/6*u_{i-1} + 2/6*u_{i}
        real, dimension(3) :: coeff1m3 = (/ 11./6.,-7./6.,  2./6./)
        ! linear weights for i - 1/2
        real, dimension(3) :: dbarcoeff = (/0.3, 0.6, 0.1/)
        

        select case (dir)
           case (1)
              imm = 2; im = 1; ip = 1; ipp = 2 
              jmm = 0; jm = 0; jp = 0; jpp = 0
           case(2)
              imm = 0; im = 0; ip = 0; ipp = 0 
              jmm = 2; jm = 1; jp = 1; jpp = 2
        end select

        do n = 1, NCONSVAR_NUMBER
           do j = jlo-1, jhi+1
               do i = ilo-1, ihi+1

                  qm2 = solnVar(n,i-imm, j-jmm)
                  qm1 = solnVar(n,i-im, j-jm)
                  q = solnVar(n,i,j)
                  qp1 = solnVar(n,i+ip, j+jp)
                  qp2 = solnVar(n,i+ipp, j+jpp)
                  
                  beta(1) = n13o12*( qm2 - 2.*qm1 + q )**2 &
                     + 0.25*( qm2 - 4.*qm1 + 3.*q )**2
                  beta(2) = n13o12*( qm1 - 2.*q + qp1 )**2 &
                     + 0.25*( qm1 - qp1 )**2
                  beta(3) = n13o12*( q - 2.*qp1 + qp2 )**2 &
                     + 0.25*( 3.*q - 4.*qp1 + qp2 )**2

                  absbeta_diff = abs(beta(3) - beta(1))

                  alpha(1:3) = dcoeff(1:3) * ( 1 + absbeta_diff/(beta(1:3) + epsilon))
                  alphabar(1:3) = dbarcoeff(1:3) * ( 1 + absbeta_diff/(beta(1:3) + epsilon)) 

                  omega(1:3) = alpha(1:3) / sum(alpha)
                  omegabar(1:3) = alphabar(1:3) / sum(alphabar)

                  qpluspoly(1) = dot_product(coeff1p1(1:3), (/qm2, qm1, q/)) 
                  qpluspoly(2) = dot_product(coeff1p2(1:3), (/qm1, q, qp1/))
                  qpluspoly(3) = dot_product(coeff1p3(1:3), (/q, qp1, qp2/))

                  qminuspoly(1) = dot_product(coeff1m1(1:3), (/qm2, qm1, q/))
                  qminuspoly(2) = dot_product(coeff1m2(1:3), (/qm1, q, qp1/))
                  qminuspoly(3) = dot_product(coeff1m3(1:3), (/q, qp1, qp2/))

                  qplus12 = dot_product(omega(1:3), qpluspoly(1:3)) 
                  qminus12 = dot_product(omegabar(1:3), qminuspoly(1:3)) 

                  fplus12 = qplus12
                  fminus12 = qminus12
                  ! monotonicity check
                  if ((qplus12 - q) * (q - qminus12) < 0.000) then
                     fplus12 = q
                     fminus12 = q
                  end if
                  if (present(check)) then
                     delp = qp1 - q
                     delm = q - qm1
                     ! tvd reconstruction
                     if (fplus12 < 0.000) then
                       fplus12 = q + 0.500 * minmod(delp, delm)
                     end if
                     if (fminus12 < 0.000) then
                       fminus12 = q - 0.500 * minmod(delp, delm)
                     end if
                  end if
                  flxminus(n,i,j) = fminus12
                  flxplus(n,i,j) = fplus12   
               end do
           end do
        end do
  end subroutine get_facevalue_weno5z

  subroutine recon_getcellfaces(dt, solnVar, x_plus, x_minus, y_plus, y_minus)

        use sim_data, only: xTpts, yTpts, dx, dy, gamma, ilo, ihi, jlo, jhi, & 
                            recon_method, iGlo, iGhi, jGlo, jGhi
        use misc_module, only: to_upper
        implicit none
        real, intent(in) :: dt
        real, dimension(1:, iGlo:, jGlo:), intent(in) :: solnVar
        real, dimension(NCONSVAR_NUMBER, iGlo:iGhi,jGlo:jGhi), intent(out) :: x_plus, x_minus, y_plus, y_minus
        
        integer :: i, j, n

        x_plus(:,:,:) = 0.0
        x_minus(:,:,:) = 0.0

        y_plus(:,:,:) = 0.0
        y_minus(:,:,:) = 0.0

        if (to_upper(trim(recon_method)) == "WENO3") then
                call get_facevalue_weno(solnVar, &
                                    x_plus, x_minus, 1)
                call get_facevalue_weno(solnVar, &
                                    y_plus, y_minus, 2)
        else if (to_upper(trim(recon_method)) == "WENO5") then
                call get_facevalue_weno5z(solnVar, &
                                    x_plus, x_minus, 1)
                call get_facevalue_weno5z(solnVar, &
                                    y_plus, y_minus, 2)
        end if

    end subroutine recon_getcellfaces

end module recon_module
