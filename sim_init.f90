module sim_init
  implicit none
contains
  subroutine init_problem(x, y, dens, velx, vely, pres, ener)
    use sim_data, only: problem, xTpts, yTpts, &
                        nx, ny, Gpts, ihi, ilo, jhi, &
                        jlo, PI, grav, usegrav
    use misc_module, only: to_upper
    use eos_module, only: eos_gete
    implicit none
    real, dimension(nx), intent(in) :: x
    real, dimension(ny), intent(in) :: y
    real, dimension(xTpts, yTpts), intent(inout) :: dens, velx, vely, pres, ener
    integer :: i, j, iInt, jInt
    
    select case (to_upper(trim(problem)))

    case ("SEDOV")
      print *, "Sedov problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           dens(i,j) = 1.0e0
           velx(i,j) = 0.0e0
           vely(i,j) = 0.0e0
           pres(i,j) = 1.0e-1
 
           ! add a region with very high pressure
           if (sqrt(x(iInt)**2 + y(jInt)**2) <= 0.1) then
             pres(i,j) = 1.0e1
           end if
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
        end do
      end do

    case ("KH")
      print *, "Kelvin-Helmholtz (KH) problem selected"
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           if (abs(y(jInt) - 0.5) > 0.25) then
             dens(i,j) = 1.0e0
             velx(i,j) = -0.5e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0
           else
             dens(i,j) = 2.0e0
             velx(i,j) = 0.5e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0
           end if
           ! give velocity perturbations
           if (abs(y(jInt) - 0.25) < 0.1) then
              vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
           if (abs(y(jInt) - 0.75) < 0.1) then
              vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
           end if
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
        end do
      end do

    case ("RT")
      print *, "Rayleigh-Taylor (RT) problem selected"
      if (.not. usegrav) then
        print *, "Gravity is needed. Setting usegrav = T"
        usegrav = .true.
      end if
      do j = jlo, jhi
        do i = ilo, ihi
           ! shift to interior index
           iInt = i - Gpts
           jInt = j - Gpts

           ! initialize the fields
           if (y(jInt) > 0.0) then
             dens(i,j) = 2.0e0
             velx(i,j) = -0.0e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0 + grav * dens(i,j) * y(jInt)
           else
             dens(i,j) = 1.0e0
             velx(i,j) = 0.0e0
             vely(i,j) = 0.0e0
             pres(i,j) = 2.5e0 + grav * dens(i,j) * y(jInt)
           end if
           ! give velocity perturbation at the interface
           !if (abs(y(jInt)) < 0.1) then
              !vely(i,j) = 0.05 * sin(2.0 * 2.0 * PI * x(iInt))
              vely(i,j) = 0.01 * (1.0 + cos(4.0 * PI * x(iInt))) &
                               * (1.0 + cos(3.0 * PI * y(jInt))) / 4.0
           !end if
           ener(i,j) = eos_gete((/dens(i,j), velx(i,j), vely(i,j), 0.0, pres(i,j)/)) 
        end do
      end do
    case default
      print *, "such problem not defined"
    end select    

     
  end subroutine init_problem
end module sim_init
