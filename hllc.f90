module hllc_module
    implicit none
contains

    subroutine hllc(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
        use eos_module, only: eos_getp
        use misc_module, only: to_upper 
        implicit none

        real(8), dimension(5), intent(in) :: Uleft, Uright
        character(len=1), intent(in) :: Dir
        real(8), dimension(5), intent(out) :: Flux
        real(8), dimension(5) :: FL, FR, UL, UR, FstarL, FstarR, UstarL, UstarR
        real(8) :: dL, udL, vdL, wdL, eL
        real(8) :: dR, udR, vdR, wdR, eR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: pL, pR, HL, HR, cL, cR
        real(8) :: ubar, vbar, wbar, Hbar, cbar
        real(8) :: qL, qR, qbar, SL, SR, SM
        real(8) :: dstarL, dstarR, pstar
        real(8) :: dustarL, dustarR, dvstarL, dvstarR
        real(8) :: dwstarL, dwstarR, estarL, estarR
        real(8), dimension(3) :: n

        dL = Uleft(1)
        udL = Uleft(2)
        vdL = Uleft(3)
        wdL = Uleft(4)
        eL = Uleft(5)

        dR = Uright(1)
        udR = Uright(2)
        vdR = Uright(3)
        wdR = Uright(4)
        eR = Uright(5)

        select case (to_upper(Dir))
        case ('X')
            n = (/1.00, 0.00, 0.00/)
        case ('Y')
            n = (/0.00, 1.00, 0.00/)
        case ('Z')
            n = (/0.00, 0.00, 1.00/)
        case default
            print *, "Wrong!! Direction should be X, Y, Z"
            stop
        end select

        sdR = sqrt(dR)
        sdL = sqrt(dL)

        ufL = udL / dL
        ufR = udR / dR

        vfL = vdL / dL
        vfR = vdR / dR

        wfL = wdL / dL
        wfR = wdR / dR

        pL =  eos_getp((/dL, ufL, vfL, wfL, eL/))
        pR =  eos_getp((/dR, ufR, vfR, wfR, eR/))

        HL = (eL + pL) / dL
        HR = (eR + pR) / dR

        cL = sqrt(gamma * pL / dL)
        cR = sqrt(gamma * pR / dR)

        ubar = (ufL * sdL + ufR * sdR) / (sdL + sdR)
        vbar = (vfL * sdL + vfR * sdR) / (sdL + sdR)
        wbar = (wfL * sdL + wfR * sdR) / (sdL + sdR)

        Hbar = (HL * sdL + HR * sdR) / (sdL + sdR)

        cbar = sqrt((gamma - 1.00) * (Hbar - 0.50 * (ubar**2 + vbar**2 + wbar**2)))

        qL = sum((/ufL, vfL, wfL/) * n)
        qR = sum((/ufR, vfR, wfR/) * n)

        qbar = sum((/ubar, vbar, wbar/) * n)

        SL = min(qL - cL, qbar - cbar)
        SR = max(qR + cR, qbar + cbar)

        SM = ((dR * qR * (SR - qR) - dL * qL * (SL - qL)) + pL - pR) / (dR * (SR - qR) - dL * (SL - qL))

        FL = (/dL * qL, dL * ufL * qL + pL * n(1), dL * vfL * qL + pL * n(2), dL * wfL * qL + pL * n(3), qL * (eL + pL)/)
        FR = (/dR * qR, dR * ufR * qR + pR * n(1), dR * vfR * qR + pR * n(2), dR * wfR * qR + pR * n(3), qR * (eR + pR)/)

        dstarL = dL * (SL - qL) / (SL - SM)
        dstarR = dR * (SR - qR) / (SR - SM)

        pstar = pL + dL * (qL - SL) * (qL - SM)

        dustarL = (dL * ufL * (SL - qL) + (pstar - pL) * n(1)) / (SL - SM)
        dustarR = (dR * ufR * (SR - qR) + (pstar - pR) * n(1)) / (SR - SM)

        dvstarL = (dL * vfL * (SL - qL) + (pstar - pL) * n(2)) / (SL - SM)
        dvstarR = (dR * vfR * (SR - qR) + (pstar - pR) * n(2)) / (SR - SM)

        dwstarL = (dL * wfL * (SL - qL) + (pstar - pL) * n(3)) / (SL - SM)
        dwstarR = (dR * wfR * (SR - qR) + (pstar - pR) * n(3)) / (SR - SM)

        estarL = ((eL * (SL - qL) - pL * qL) + pstar * SM) / (SL - SM)
        estarR = ((eR * (SR - qR) - pR * qR) + pstar * SM) / (SR - SM)

        UstarL = (/dstarL, dustarL, dvstarL, dwstarL, estarL/)
        UstarR = (/dstarR, dustarR, dvstarR, dwstarR, estarR/)

        UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL/)
        UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR/)

        FstarL = FL + SL * (UstarL - UL)
        FstarR = FR + SR * (UstarR - UR)

        if (0.000 < SL) then
            Flux = FL
        else if (SL <= 0.000 .and. 0.000 < SM) then
            Flux = FstarL
        else if (SM <= 0.000 .and. 0.000 <= SR) then
            Flux = FstarR
        else
            Flux = FR
        end if
        return
    end subroutine hllc

end module hllc_module

