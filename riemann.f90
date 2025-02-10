module riemann_module
#include "header.h"
    implicit none

    ! HLLC solver is only for pure hydro simulation
    ! HLLE can be used for both pure Hydro and MHD simulation

contains

#ifndef MHD
    subroutine hllc(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
        use eos_module, only: eos_gete
        use misc_module, only: to_upper 
        implicit none

        real(8), dimension(5), intent(in) :: Uleft, Uright
        character(len=1), intent(in) :: Dir
        real(8), dimension(5), intent(out) :: Flux
        real(8), dimension(5) :: FL, FR, UL, UR, FstarL, FstarR, UstarL, UstarR
        real(8) :: dL, udL, vdL, wdL, pL
        real(8) :: dR, udR, vdR, wdR, pR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, HL, HR, cL, cR
        real(8) :: ubar, vbar, wbar, Hbar, cbar
        real(8) :: qL, qR, qbar, SL, SR, SM
        real(8) :: dstarL, dstarR, pstar
        real(8) :: dustarL, dustarR, dvstarL, dvstarR
        real(8) :: dwstarL, dwstarR, estarL, estarR
        real(8), dimension(3) :: n

        dL = Uleft(1)
        ufL = Uleft(2)
        vfL = Uleft(3)
        wfL = Uleft(4)
        pL = Uleft(5)

        dR = Uright(1)
        ufR = Uright(2)
        vfR = Uright(3)
        wfR = Uright(4)
        pR = Uright(5)

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

        udL = ufL * dL
        udR = ufR * dR

        vdL = vfL * dL
        vdR = vfR * dR

        wdL = wfL * dL
        wdR = wfR * dR

        eL =  eos_gete(Uleft(1:5))
        eR =  eos_gete(Uright(1:5))

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

        FL = (/dL * qL, &
               dL * ufL * qL + pL * n(1), &
               dL * vfL * qL + pL * n(2), &
               dL * wfL * qL + pL * n(3), &
               qL * (eL + pL)/)
        FR = (/dR * qR, &
               dR * ufR * qR + pR * n(1), &
               dR * vfR * qR + pR * n(2), &
               dR * wfR * qR + pR * n(3), &
               qR * (eR + pR)/)

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
#endif

    subroutine hlle(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
#ifdef MHD
        use sim_data, only: ch
#endif
        use eos_module, only: eos_gete
        use misc_module, only: to_upper 
        implicit none

        character(len=1), intent(in) :: Dir
#ifdef MHD
        real(8), dimension(9), intent(in) :: Uleft, Uright
        real(8), dimension(9), intent(out) :: Flux
        real(8), dimension(9) :: FL, FR, UL, UR, Fstar
#else
        real(8), dimension(5), intent(in) :: Uleft, Uright
        real(8), dimension(5), intent(out) :: Flux
        real(8), dimension(5) :: FL, FR, UL, UR, Fstar
#endif
        real(8) :: dL, udL, vdL, wdL, pL, bxL, byL, bzL, bpL
        real(8) :: dR, udR, vdR, wdR, pR, bxR, byR, bzR, bpR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, HL, HR, cL, cR, B2L, B2R, vBL, vBR
        real(8) :: ubar, vbar, wbar, Hbar, cbar
        real(8) :: qL, qR, qbar, qBL, qBR, SL, SR, cfL, cfR
        real(8), dimension(3) :: n

        dL = Uleft(1)
        ufL = Uleft(2)
        vfL = Uleft(3)
        wfL = Uleft(4)
        pL = Uleft(5)

        dR = Uright(1)
        ufR = Uright(2)
        vfR = Uright(3)
        wfR = Uright(4)
        pR = Uright(5)

#ifdef MHD
        bxL = Uleft(6)
        byL = Uleft(7)
        bzL = Uleft(8)
        bpL = Uleft(9)
        
        bxR = Uright(6)
        byR = Uright(7)
        bzR = Uright(8)
        bpR = Uright(9)
#endif

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

        udL = ufL * dL
        udR = ufR * dR

        vdL = vfL * dL
        vdR = vfR * dR

        wdL = wfL * dL
        wdR = wfR * dR

#ifdef MHD
        eL =  eos_gete(Uleft(1:8))
        eR =  eos_gete(Uright(1:8))
#else
        eL =  eos_gete(Uleft(1:5))
        eR =  eos_gete(Uright(1:5))
#endif

        cL = sqrt(gamma * pL / dL)
        cR = sqrt(gamma * pR / dR)

        qL = sum((/ufL, vfL, wfL/) * n)
        qR = sum((/ufR, vfR, wfR/) * n)

#ifndef MHD
        HL = (eL + pL) / dL
        HR = (eR + pR) / dR

        ubar = (ufL * sdL + ufR * sdR) / (sdL + sdR)
        vbar = (vfL * sdL + vfR * sdR) / (sdL + sdR)
        wbar = (wfL * sdL + wfR * sdR) / (sdL + sdR)

        Hbar = (HL * sdL + HR * sdR) / (sdL + sdR)

        cbar = sqrt((gamma - 1.00) * (Hbar - 0.50 * (ubar**2 + vbar**2 + wbar**2)))

        qbar = sum((/ubar, vbar, wbar/) * n)

        SL = min(min(qL - cL, qbar - cbar), 0.0)
        SR = max(max(qR + cR, qbar + cbar), 0.0)

        FL = (/dL * qL, &
               dL * ufL * qL + pL * n(1), &
               dL * vfL * qL + pL * n(2), &
               dL * wfL * qL + pL * n(3), &
               qL * (eL + pL)/)
        FR = (/dR * qR, &
               dR * ufR * qR + pR * n(1), &
               dR * vfR * qR + pR * n(2), &
               dR * wfR * qR + pR * n(3), &
               qR * (eR + pR)/)

        UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL/)
        UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR/)

        Fstar = (SR * FL - SL * FR + (SL * SR * (UR - UL))) / (SR - SL)
#endif
#ifdef MHD
        ! magnetic field value in the normal direction
        qBL = sum((/bxL, byL, bzL/) * n)
        qBR = sum((/bxR, byR, bzR/) * n)
        
        ! B.B quantity
        B2L = sum(Uleft(6:8) * Uleft(6:8))
        B2R = sum(Uright(6:8) * Uright(6:8))

        ! V.B quantity
        vBL = sum(Uleft(2:4) * Uleft(6:8))
        vBR = sum(Uright(2:4) * Uright(6:8))
        
        ! fast magnetosonic wave speed
        ! Miyoshi eqn 67 HLLD paper 
        cfL = sqrt(((gamma * pL + B2L) + &
          sqrt((gamma * pL + B2L)**2 - gamma * pL * 4.0 * qBL**2))/(2.0 * dL))

        cfR = sqrt(((gamma * pR + B2R) + &
          sqrt((gamma * pR + B2R)**2 - gamma * pR * 4.0 * qBR**2))/(2.0 * dR))

        SL = min(qL, qR) - max(cfL, cfR) 
        SR = max(qL, qR) + max(cfL, cfR) 
        
        FL = (/dL * qL, &
               dL * ufL * qL - qBL * bxL + (0.50 * B2L + pL) * n(1) , &
               dL * vfL * qL - qBL * byL + (0.50 * B2L + pL) * n(2) , &
               dL * wfL * qL - qBL * bzL + (0.50 * B2L + pL) * n(3) , &
               qL * (eL + pL + 0.50 * B2L) - qBL * vBL              , &
               qL * bxL - qBL * ufL                                 , &
               qL * byL - qBL * vfL                                 , &
               qL * bzL - qBL * wfL                                 , &
               0.0/)
        FR = (/dR * qR, &
               dR * ufR * qR - qBR * bxR + (0.50 * B2R + pR) * n(1) , &
               dR * vfR * qR - qBR * byR + (0.50 * B2R + pR) * n(2) , &
               dR * wfR * qR - qBR * bzR + (0.50 * B2R + pR) * n(3) , &
               qR * (eR + pR + 0.50 * B2R) - qBR * vBR              , &
               qR * bxR - qBR * ufR                                 , &
               qR * byR - qBR * vfR                                 , &
               qR * bzR - qBR * wfR                                 , &
               0.0/)

        UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL, bxL, byL, bzL, bpL/)
        UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR, bxR, byR, bzR, bpR/)

        Fstar = (SR * FL - SL * FR + (SL * SR * (UR - UL))) / (SR - SL)
#endif
        if (0.000 < SL) then
            Flux = FL
        else if (SR < 0.000) then
            Flux = FR
        else
            Flux = Fstar
        end if
#ifdef MHD
        Flux(9) = 0.50 * ((ch * ch) * (qBL + qBR) - ch * (bpR - bpL))
        select case (to_upper(Dir))
        case ('X')
           Flux(6) = 0.50 * ((bpL + bpR) - ch * (qBR - qBL))
        case ('Y')
           Flux(7) = 0.50 * ((bpL + bpR) - ch * (qBR - qBL))
        case ('Z')
           Flux(8) = 0.50 * ((bpL + bpR) - ch * (qBR - qBL))
        end select
#endif
        return
    end subroutine hlle

#ifdef MHD
    subroutine hlld(Uleft, Uright, Dir, Flux)
        use sim_data, only: gamma
        use eos_module, only: eos_gete
        use misc_module, only: to_upper 
        implicit none

        real(8), dimension(9), intent(in) :: Uleft, Uright
        character(len=1), intent(in) :: Dir
        real(8), dimension(9), intent(out) :: Flux
        real(8), dimension(9) :: FL, FR, UL, UR, FstarL, FstarR, UstarL, UstarR
        real(8) :: dL, udL, vdL, wdL, pL, bxL, byL, bzL, bpL
        real(8) :: dR, udR, vdR, wdR, pR, bxR, byR, bzR, bpR
        real(8) :: sdR, sdL, ufL, ufR, vfL, vfR, wfL, wfR
        real(8) :: eL, eR, cfL, cfR, ptL, ptR 
        real(8) :: qL, qR, SL, SR, SM, qBL, qBR, B2L, B2R, vBL, vBR
        real(8) :: dstarL, dstarR, pstar
        real(8) :: dustarL, dustarR, dvstarL, dvstarR
        real(8) :: dwstarL, dwstarR, estarL, estarR
        real(8), dimension(3) :: n

        dL = Uleft(1)
        ufL = Uleft(2)
        vfL = Uleft(3)
        wfL = Uleft(4)
        pL = Uleft(5)
        bxL = Uleft(6)
        byL = Uleft(7)
        bzL = Uleft(8)
        bpL = Uleft(9)

        dR = Uright(1)
        ufR = Uright(2)
        vfR = Uright(3)
        wfR = Uright(4)
        pR = Uright(5)
        bxR = Uright(6)
        byR = Uright(7)
        bzR = Uright(8)
        bpR = Uright(9)

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

        udL = ufL * dL
        udR = ufR * dR

        vdL = vfL * dL
        vdR = vfR * dR

        wdL = wfL * dL
        wdR = wfR * dR

        eL =  eos_gete(Uleft(1:8))
        eR =  eos_gete(Uright(1:8))

        qL = sum((/ufL, vfL, wfL/) * n)
        qR = sum((/ufR, vfR, wfR/) * n)

        ! magnetic field value in the normal direction
        qBL = sum((/bxL, byL, bzL/) * n)
        qBR = sum((/bxR, byR, bzR/) * n)
        
        ! B.B quantity
        B2L = sum(Uleft(6:8) * Uleft(6:8))
        B2R = sum(Uright(6:8) * Uright(6:8))

        ! V.B quantity
        vBL = sum(Uleft(2:4) * Uleft(6:8))
        vBR = sum(Uright(2:4) * Uright(6:8))

        ! fast magnetosonic wave speed
        ! Miyoshi eqn 67 HLLD paper 
        cfL = sqrt(((gamma * pL + B2L) + &
          sqrt((gamma * pL + B2L)**2 - gamma * pL * 4.0 * qBL**2))/(2.0 * dL))

        cfR = sqrt(((gamma * pR + B2R) + &
          sqrt((gamma * pR + B2R)**2 - gamma * pR * 4.0 * qBR**2))/(2.0 * dR))

        ptL = pL + 0.5e0 * B2L
        ptR = pR + 0.5e0 * B2R

        FL = (/dL * qL, &
               dL * ufL * qL - qBL * bxL + (0.50 * B2L + pL) * n(1) , &
               dL * vfL * qL - qBL * byL + (0.50 * B2L + pL) * n(2) , &
               dL * wfL * qL - qBL * bzL + (0.50 * B2L + pL) * n(3) , &
               qL * (eL + pL + 0.50 * B2L) - qBL * vBL              , &
               qL * bxL - qBL * ufL                                 , &
               qL * byL - qBL * vfL                                 , &
               qL * bzL - qBL * wfL                                 , &
               0.0/)
        FR = (/dR * qR, &
               dR * ufR * qR - qBR * bxR + (0.50 * B2R + pR) * n(1) , &
               dR * vfR * qR - qBR * byR + (0.50 * B2R + pR) * n(2) , &
               dR * wfR * qR - qBR * bzR + (0.50 * B2R + pR) * n(3) , &
               qR * (eR + pR + 0.50 * B2R) - qBR * vBR              , &
               qR * bxR - qBR * ufR                                 , &
               qR * byR - qBR * vfR                                 , &
               qR * bzR - qBR * wfR                                 , &
               0.0/)

        SL = min(qL, qR) - max(cfL, cfR) 
        SR = max(qL, qR) + max(cfL, cfR) 
        
        UL = (/dL, dL * ufL, dL * vfL, dL * wfL, eL, bxL, byL, bzL, bpL/)
        UR = (/dR, dR * ufR, dR * vfR, dR * wfR, eR, bxR, byR, bzR, bpR/)
        
        if (0.000 > SL) then
            Flux = FL
        else if (SR < 0.000) then
            Flux = FR
        else
           ! to do
           Flux = 0.0 
        end if
        return
    end subroutine hlld
#endif


end module riemann_module

