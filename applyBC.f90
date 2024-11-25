module applyBC_module
  implicit none
contains

  subroutine applyBC(q, Dir, flip)
    use misc_module, only: to_upper
    use sim_data, only: xbctype, ybctype, &
                        xTpts, yTpts
    implicit none
    real(8), dimension(xTpts, yTpts), intent(inout) :: q
    character(len=1), intent(in) :: Dir
    logical, optional, intent(in) :: flip
    real :: sig
   
    select case (to_upper(Dir))
    case ('X')
      if (to_upper(trim(xbctype)) == "PERIODIC") then
        q(1, :) = q(xTpts-3, :)
        q(2, :) = q(xTpts-2, :)
        q(xTpts-1, :) = q(3, :)
        q(xTpts, :) = q(4, :)
      else if (to_upper(trim(xbctype)) == "REFLECT") then
        sig = 1.000
        if (present(flip)) sig = -1.000
        q(1, :) = sig * q(4, :)
        q(2, :) = sig * q(3, :)
        q(xTpts-1, :) = sig * q(xTpts-2, :)
        q(xTpts, :) = sig * q(xTpts-3, :)
      end if 
    case ('Y')
      if (to_upper(trim(ybctype)) == "PERIODIC") then
        q(:, 1) = q(:, yTpts-3)
        q(:, 2) = q(:, yTpts-2)
        q(:, yTpts-1) = q(:, 3)
        q(:, yTpts) = q(:, 4)
      else if (to_upper(trim(ybctype)) == "REFLECT") then
        sig = 1.000
        if (present(flip)) sig = -1.000
        q(:, 1) = sig * q(:, 4)
        q(:, 2) = sig * q(:, 3)
        q(:, yTpts-1) = sig * q(:, yTpts-2)
        q(:, yTpts) = sig * q(:, yTpts-3)
      end if 
    case default
        print *, "Wrong!! This code only solves in 2D"
        stop
    end select
  end subroutine applyBC
  
  subroutine applyBC_all(dens, velx, vely, pres, ener)
    use sim_data, only: xTpts, yTpts
    implicit none
    real, dimension(xTpts, yTpts), intent(inout) :: dens, velx, vely, pres, ener

    call applyBC(dens, "x") 
    call applyBC(velx, "x", .true.) 
    call applyBC(vely, "x") 
    call applyBC(pres, "x") 
    call applyBC(ener, "x") 
    
    call applyBC(dens, "y") 
    call applyBC(velx, "y") 
    call applyBC(vely, "y", .true.) 
    call applyBC(pres, "y") 
    call applyBC(ener, "y")

  end subroutine applyBC_all
end module applyBC_module

