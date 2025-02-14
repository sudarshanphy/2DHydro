module sim_restart
#include "param.h"
  implicit none
contains
  subroutine restart_problem(fno, time, solnVar)
    use sim_data, only: xTpts, yTpts, basenm, ilo, ihi, jlo, jhi
    implicit none
    integer, intent(in) :: fno
    real, intent(out) :: time
    real, dimension(xTpts, yTpts, NVAR_NUMBER), intent(out) :: solnVar
    integer :: i, j, ionum, pos1, pos2
    character(len=256) :: fname, line1, line2, line3, line4, line5, linehash
    character(len=10) :: int_to_str 
    real :: x, y
    
    solnVar(:,:,:) = 0.0

    ionum = 56 
12 format (1x, 50(es25.18, :, 1x))

   write(int_to_str, "(I4.4)") fno
   fname = "./output/"//trim(adjustl(basenm))//"_"//trim(adjustl(int_to_str))//".dat"
   open(unit=ionum, file=fname, status="old")
   read(ionum, '(A)') line1
   read(ionum, *) line2
   read(ionum, *) line3
   read(ionum, *) line4
   read(ionum, *) line5
    
   do i = ilo, ihi
     do j = jlo, jhi
         read(ionum, 12) x, y, solnVar(i,j,NVAR_BEGIN:NVAR_END)
     end do
     read(ionum, *) linehash
   end do 
   close(ionum)

   ! Find the first instance of whitespace.  Split label and data.
   pos1 = scan(line1, '=')
   pos2 = scan(line1, ",")
   read(line1(pos1+1:pos2-1), *) time

  end subroutine restart_problem
end module sim_restart
