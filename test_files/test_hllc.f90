program test_hllc
    use riemann_module, only: hllc

    implicit none
    real, dimension(5) :: Uleft, Uright, F
    real :: gamma
    character(len=1) :: Dir

    Uleft = (/1.0,0.0,0.0,0.0,2.5/)
    Uright = (/1.0,0.0,0.0,0.0,2.5/)

    gamma = 1.4
    Dir = "z"

    call hllc(Uleft, Uright, gamma, Dir, F)
    print *, "F = ", F

end program test_hllc
