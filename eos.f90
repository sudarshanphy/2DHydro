module eos_module
  implicit none
contains

    function eos_getp(U) result(pf)
        use sim_data, only: gamma
        implicit none
        real :: U(5)
        real :: pf
        real :: df, uf, vf, wf, ef

        df = U(1)
        uf = U(2)
        vf = U(3)
        wf = U(4)
        ef = U(5)

        pf = (ef - 0.5 * df * (uf**2 + vf**2 + wf**2)) * (gamma - 1.0)
    end function eos_getp


    function eos_gete(U) result(ef)
        use sim_data, only: gamma
        implicit none
        real :: U(5)
        real :: ef
        real :: df, uf, vf, wf, pf

        df = U(1)
        uf = U(2)
        vf = U(3)
        wf = U(4)
        pf = U(5)

        ef = pf/(gamma - 1.0) + 0.5 * df * (uf**2 + vf**2 + wf**2)
    end function eos_gete

end module eos_module
