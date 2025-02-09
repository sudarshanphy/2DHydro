module recon_module
    implicit none
contains
  function minmod(a,b) result(c)
    implicit none
    real :: a, b
    real :: c, fact

    fact = sign(1.0, a) + sign(1.0, b)
    c = 0.500 * fact * min(abs(a), abs(b)) 

  end function minmod

  subroutine get_facevalue_weno(qm1, q, qp1, fplus12, fminus12, check)
        implicit none
        real, intent(in) :: qm1, q, qp1
        real, intent(out):: fplus12, fminus12
        logical, intent(in), optional :: check

        real, parameter :: epsilon = 1.0e-6
        real, dimension(2) :: beta, alpha, alphabar, &
                                 omega, omegabar, & 
                                 qpluspoly, qminuspoly
        real :: qplus12, qminus12
        integer :: i
        real :: delp, delm
        ! coefficients for the polynomial
        real, dimension(2) :: c0 = (/1.00/2.00, 1.00/2.00/)
        real, dimension(2) :: c1 = (/-1.00/2.00, 3.00/2.00/)
        real, dimension(2) :: c_1 = (/3.00/2.00, -1.00/2.00/)
        ! linear weights 
        real, dimension(2) :: dcoeff = (/2.00/3.00, 1.00/3.00/)
        real, dimension(2) :: dbarcoeff = (/1.00/3.00, 2.00/3.00/)

        beta(1) = (qp1 - q)**2
        beta(2) = (q - qm1)**2

        do i = 1, 2
            alpha(i) = dcoeff(i) / (beta(i) + epsilon)**2
            alphabar(i) = dbarcoeff(i) / (beta(i) + epsilon)**2
        end do

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

  end subroutine get_facevalue_weno

  subroutine get_facevalue_weno5z(qm2, qm1, q, qp1, qp2, fplus12, fminus12, check)
        implicit none
        real, intent(in) :: qm2, qm1, q, qp1, qp2
        real, intent(out):: fplus12, fminus12
        logical, intent(in), optional :: check

        real, parameter :: epsilon = 1.0e-36
        real, dimension(3) :: beta, alpha, alphabar, &
                              omega, omegabar, qpluspoly, &
                              qminuspoly
        real :: qplus12, qminus12
        integer :: i
        real :: delp, delm, absbeta_diff

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
      

        beta(1) = n13o12*( qm2 - 2.*qm1 + q )**2 &
           + 0.25*( qm2 - 4.*qm1 + 3.*q )**2
        beta(2) = n13o12*( qm1 - 2.*q + qp1 )**2 &
           + 0.25*( qm1 - qp1 )**2
        beta(3) = n13o12*( q - 2.*qp1 + qp2 )**2 &
           + 0.25*( 3.*q - 4.*qp1 + qp2 )**2

         absbeta_diff = abs(beta(3) - beta(1))

        do i = 1, 3
            alpha(i) = dcoeff(i) * ( 1 + absbeta_diff/(beta(i) + epsilon))
            alphabar(i) = dbarcoeff(i) * ( 1 + absbeta_diff/(beta(i) + epsilon)) 
        end do

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

  end subroutine get_facevalue_weno5z

  subroutine recon_getcellfaces(dens, velx, vely, pres, dt, &
                                   !method, slimiter, &
                                   xr_plus, xu_plus, xv_plus, xp_plus, &
                                   xr_minus, xu_minus, xv_minus, xp_minus, &
                                   yr_plus, yu_plus, yv_plus, yp_plus, &
                                   yr_minus, yu_minus, yv_minus, yp_minus)

        use sim_data, only: xTpts, yTpts, dx, dy, gamma, ilo, ihi, jlo, jhi, & 
                            recon_method
        use misc_module, only: to_upper
        implicit none
        real, intent(in) :: dt
        real, dimension(xTpts,yTpts), intent(in) :: dens, velx, vely, pres
        !character(len=*) :: method
        !character(len=*) :: slimiter

        real, dimension(xTpts, yTpts), intent(out) :: xr_plus, xu_plus, xv_plus, xp_plus
        real, dimension(xTpts, yTpts), intent(out) :: xr_minus, xu_minus, xv_minus, xp_minus
        real, dimension(xTpts, yTpts), intent(out) :: yr_plus, yu_plus, yv_plus, yp_plus
        real, dimension(xTpts, yTpts), intent(out) :: yr_minus, yu_minus, yv_minus, yp_minus

        integer :: i, j

        xr_plus = 0.0
        xu_plus = 0.0
        xv_plus = 0.0
        xp_plus = 0.0

        xr_minus = 0.0
        xu_minus = 0.0
        xv_minus = 0.0
        xp_minus = 0.0

        yr_plus = 0.0
        yu_plus = 0.0
        yv_plus = 0.0
        yp_plus = 0.0

        yr_minus = 0.0
        yu_minus = 0.0
        yv_minus = 0.0
        yp_minus = 0.0

        if (to_upper(trim(recon_method)) == "WENO3") then
            do i = ilo-1, ihi+1
                do j = jlo-1, jhi+1

                    call get_facevalue_weno(dens(i - 1, j), dens(i, j), dens(i + 1, j), &
                      xr_plus(i, j), xr_minus(i, j))
                    call get_facevalue_weno(velx(i - 1, j), velx(i, j), velx(i + 1, j), &
                      xu_plus(i, j), xu_minus(i, j))
                    call get_facevalue_weno(vely(i - 1, j), vely(i, j), vely(i + 1, j), &
                      xv_plus(i, j), xv_minus(i, j))
                    call get_facevalue_weno(pres(i - 1, j), pres(i, j), pres(i + 1, j), &
                      xp_plus(i, j), xp_minus(i, j))

                    call get_facevalue_weno(dens(i, j - 1), dens(i, j), dens(i, j + 1), &
                      yr_plus(i, j), yr_minus(i, j))
                    call get_facevalue_weno(velx(i, j - 1), velx(i, j), velx(i, j + 1), &
                      yu_plus(i, j), yu_minus(i, j))
                    call get_facevalue_weno(vely(i, j - 1), vely(i, j), vely(i, j + 1), &
                      yv_plus(i, j), yv_minus(i, j))
                    call get_facevalue_weno(pres(i, j - 1), pres(i, j), pres(i, j + 1), &
                      yp_plus(i, j), yp_minus(i, j))

                end do
            end do
        else if (to_upper(trim(recon_method)) == "WENO5") then
            do i = ilo-1, ihi+1
                do j = jlo-1, jhi+1

                call get_facevalue_weno5z(dens(i - 2, j), dens(i - 1, j), dens(i, j), &
                  dens(i + 1, j), dens(i + 2, j), xr_plus(i, j), xr_minus(i, j))
                call get_facevalue_weno5z(velx(i - 2, j), velx(i - 1, j), velx(i, j), &
                  velx(i + 1, j), velx(i + 2, j), xu_plus(i, j), xu_minus(i, j))
                call get_facevalue_weno5z(vely(i - 2, j), vely(i - 1, j), vely(i, j), &
                  vely(i + 1, j), vely(i + 2, j), xv_plus(i, j), xv_minus(i, j))
                call get_facevalue_weno5z(pres(i - 2, j), pres(i - 1, j), pres(i, j), &
                  pres(i + 1, j), pres(i + 2, j), xp_plus(i, j), xp_minus(i, j))
                                                                                                                   
                call get_facevalue_weno5z(dens(i, j - 2), dens(i, j - 1), dens(i, j), &
                  dens(i, j + 1), dens(i, j + 2), yr_plus(i, j), yr_minus(i, j))
                call get_facevalue_weno5z(velx(i, j - 2), velx(i, j - 1), velx(i, j), &
                  velx(i, j + 1), velx(i, j + 2), yu_plus(i, j), yu_minus(i, j))
                call get_facevalue_weno5z(vely(i, j - 2), vely(i, j - 1), vely(i, j), &
                  vely(i, j + 1), vely(i, j + 2), yv_plus(i, j), yv_minus(i, j))
                call get_facevalue_weno5z(pres(i, j - 2), pres(i, j - 1), pres(i, j), &
                  pres(i, j + 1), pres(i, j + 2), yp_plus(i, j), yp_minus(i, j))

                end do
            end do
        end if

    end subroutine recon_getcellfaces

end module recon_module
