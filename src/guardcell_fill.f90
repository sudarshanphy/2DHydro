subroutine guardcell_fill()
#include "param.h"
  use sim_data, only: ilo, jlo, ihi, jhi, &
                      iGlo, jGlo, iGhi, jGhi, &
                      mainVar, myrank, xblk, yblk, &
                      lxTpts, lyTpts, comm, ierr, status1, &
                      at_xlboundary, at_xrboundary, at_ylboundary, &
                      at_yrboundary, mainVar, Gpts, dtype_mpi

  use mpi
  implicit none
  integer :: l_blk, r_blk, t_blk, b_blk
  integer :: tag1 = 20, tag2 = 21, tag3 = 23, tag4 = 24
  integer :: xnum, ynum
  integer :: i, j, k, count


  ! see if we don't have neighbouring blocks
  l_blk = merge(MPI_PROC_NULL, myrank - 1, at_xlboundary)
  r_blk = merge(MPI_PROC_NULL, myrank + 1, at_xrboundary)
  t_blk = merge(MPI_PROC_NULL, myrank + xblk, at_yrboundary)
  b_blk = merge(MPI_PROC_NULL, myrank - xblk, at_ylboundary)
 
  if (xblk > 1) then
    ! communication in X-direction
    
    ! send data to right
    call MPI_SENDRECV(mainVar(1,ihi-Gpts+1,jGlo), &
                      1, dtype_mpi(1), r_blk, tag1, &
                      mainVar(1,iGlo,jGlo), &
                      1, dtype_mpi(1), l_blk, tag1, &
                      comm, status1, ierr)

    ! send data to left
    call MPI_SENDRECV(mainVar(1,ilo,jGlo), &
                      1, dtype_mpi(1), l_blk, tag2, &
                      mainVar(1,ihi+1,jGlo), &
                      1, dtype_mpi(1), r_blk, tag2, &
                      comm, status1, ierr)
    
  endif
  if (yblk > 1) then
    ! communication in Y-direction
    
    ! send data to top
    call MPI_SENDRECV(mainVar(1,iGlo,jhi-Gpts+1), &
                      1, dtype_mpi(2), t_blk, tag3, &
                      mainVar(1,iGlo,jGlo), &
                      1, dtype_mpi(2), b_blk, tag3, &
                      comm, status1, ierr)
    
    ! send data to bottom
    call MPI_SENDRECV(mainVar(1,iGlo,jlo), &
                      1, dtype_mpi(2), b_blk, tag4, &
                      mainVar(1,iGlo,jhi+1), &
                      1, dtype_mpi(2), t_blk, tag4, &
                      comm, status1, ierr)
  endif
end subroutine guardcell_fill
