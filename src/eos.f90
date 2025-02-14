module eos_module
#include "header.h"
#include "param.h"
  implicit none
contains

    function eos_getp(U) result(pf)
        use sim_data, only: gamma
        implicit none
        real, dimension(NVAR_NUMBER) :: U
        real :: pf
        real :: df, uf, vf, wf, ef, bxf, byf, bzf

        df = U(DENS_VAR)
        uf = U(VELX_VAR)
        vf = U(VELY_VAR)
        wf = U(VELZ_VAR)
        ef = U(ENER_VAR)
        bxf = 0.0; byf = 0.0; bzf = 0.0
#ifdef MHD
        bxf = U(BMFX_VAR)
        byf = U(BMFY_VAR)
        bzf = U(BMFZ_VAR)
#endif

        pf = (ef - 0.5 * (df * (uf**2 + vf**2 + wf**2) + &
                               (bxf**2 + byf**2 + bzf**2))) * (gamma - 1.0)
    end function eos_getp


    function eos_gete(U) result(ef)
        use sim_data, only: gamma
        implicit none
        real, dimension(NVAR_NUMBER) :: U
        real :: ef
        real :: df, uf, vf, wf, pf, bxf, byf, bzf

        df = U(DENS_VAR)
        uf = U(VELX_VAR)
        vf = U(VELY_VAR)
        wf = U(VELZ_VAR)
        pf = U(PRES_VAR)
        bxf = 0.0; byf = 0.0; bzf = 0.0
#ifdef MHD
        bxf = U(BMFX_VAR)
        byf = U(BMFY_VAR)
        bzf = U(BMFZ_VAR)
#endif

        ef = pf/(gamma - 1.0) + 0.5 * (df * (uf**2 + vf**2 + wf**2) + &
                                         (bxf**2 + byf**2 + bzf**2))
    end function eos_gete

end module eos_module
