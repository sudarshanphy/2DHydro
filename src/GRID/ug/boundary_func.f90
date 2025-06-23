module boundary_func
#include "param.h"

  implicit none
  contains

  subroutine boundary_init()
    use sim_data
    use misc_module, only: to_upper
    implicit none

    ! is the block at the boundary
    at_xlboundary = .false.
    at_xrboundary = .false.
    at_ylboundary = .false.
    at_yrboundary = .false.

    if (ilo == 1)  at_xlboundary = .true.
    if (ihi == nx) at_xrboundary = .true.
    if (jlo == 1)  at_ylboundary = .true.
    if (jhi == ny) at_yrboundary = .true.

    if (to_upper(trim(xlbctype)) == "FLOW") then 
        xlbc_int = FLOW
    else if (to_upper(trim(xlbctype)) == "REFLECT") then 
        xlbc_int = REFLECT
    else if (to_upper(trim(xlbctype)) == "PERIODIC") then
        xlbc_int = PERIODIC
    else if (to_upper(trim(xlbctype)) == "CUSTOM") then
        xlbc_int = CUSTOM
    else
        print *, "X-left BC options are: flow, reflect, periodic, custom "
        stop
    end if

     if (to_upper(trim(xrbctype)) == "FLOW") then
        xrbc_int = FLOW
    else if (to_upper(trim(xrbctype)) == "REFLECT") then 
        xrbc_int = REFLECT
    else if (to_upper(trim(xrbctype)) == "PERIODIC") then
        xrbc_int = PERIODIC
    else if (to_upper(trim(xrbctype)) == "CUSTOM") then
        xrbc_int = CUSTOM
    else
        print *, "X-right BC options are: flow, reflect, periodic, custom "
        stop
    end if

    if (to_upper(trim(ylbctype)) == "FLOW") then
        ylbc_int = FLOW
    else if (to_upper(trim(ylbctype)) == "REFLECT") then
        ylbc_int = REFLECT
    else if (to_upper(trim(ylbctype)) == "PERIODIC") then
        ylbc_int = PERIODIC
    else if (to_upper(trim(ylbctype)) == "CUSTOM") then
        ylbc_int = CUSTOM
    else
        print *, "Y-left BC options are: flow, reflect, periodic, custom "
        stop
    end if

     if (to_upper(trim(yrbctype)) == "FLOW") then
        yrbc_int = FLOW
    else if (to_upper(trim(yrbctype)) == "REFLECT") then 
        yrbc_int = REFLECT
    else if (to_upper(trim(yrbctype)) == "PERIODIC") then
        yrbc_int = PERIODIC
    else if (to_upper(trim(yrbctype)) == "CUSTOM") then
        yrbc_int = CUSTOM
    else
        print *, "Y-right BC options are: flow, reflect, periodic, custom "
        stop
    end if
 end subroutine boundary_init

end module boundary_func
