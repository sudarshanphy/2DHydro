module eos_module
#include "header.h"
  implicit none
contains

    function eos_getp(U) result(pf)
        use sim_data, only: gamma
        implicit none
#ifdef MHD
        real :: U(9)
#else
        real :: U(5)
#endif
        real :: pf
        real :: df, uf, vf, wf, ef, bxf, byf, bzf

        df = U(1)
        uf = U(2)
        vf = U(3)
        wf = U(4)
        ef = U(5)
        bxf = 0.0; byf = 0.0; bzf = 0.0
#ifdef MHD
        bxf = U(6)
        byf = U(7)
        bzf = U(8)
#endif

        pf = (ef - 0.5 * (df * (uf**2 + vf**2 + wf**2) + &
                               (bxf**2 + byf**2 + bzf**2))) * (gamma - 1.0)
    end function eos_getp


    function eos_gete(U) result(ef)
        use sim_data, only: gamma
        implicit none
#ifdef MHD
        real :: U(9)
#else
        real :: U(5)
#endif
        real :: ef
        real :: df, uf, vf, wf, pf, bxf, byf, bzf

        df = U(1)
        uf = U(2)
        vf = U(3)
        wf = U(4)
        pf = U(5)
        bxf = 0.0; byf = 0.0; bzf = 0.0
#ifdef MHD
        bxf = U(6)
        byf = U(7)
        bzf = U(8)
#endif

        ef = pf/(gamma - 1.0) + 0.5 * (df * (uf**2 + vf**2 + wf**2) + &
                                         (bxf**2 + byf**2 + bzf**2))
    end function eos_gete

end module eos_module
