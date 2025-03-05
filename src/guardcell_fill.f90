subroutine guardcell_fill(solnVar)
#include "param.h"
  use sim_data, only: ilo, jlo, ihi, jhi, &
                      iGlo, jGlo, iGhi, jGhi, &
                      mainVar, myrank, xblk, yblk, &
                      lxTpts, lyTpts, comm, ierr, status1

  use mpi
  implicit none
  real, pointer :: solnVar(:,:,:)
  integer :: l_blk, r_blk, t_blk, b_blk
  integer :: sendtag = 20, recvtag = 21

  l_blk = myrank - 1
  r_blk = myrank + 1
  t_blk = myrank + xblk
  b_blk = myrank - xblk

  if (mod(l_blk,xblk) >= 0) then
    call MPI_SEND(solnVar(ilo:ilo+Gpts-1,:,:), Gpts*lyTpts*NVAR_NUMBER, MPI_DOUBLE, &
                  l_blk, sendtag, comm, ierr)
    call MPI_RECV(solnVar(iGlo:ilo-1,:,:), Gpts*lyTpys*NVAR_NUMBER, MPI_DOUBLE, &
                  l_blk, recvtag, comm, status1, ierr)
  endif
  if (mod(r_blk,xblk) <= xblk) then
    call MPI_SEND(solnVar(ihi-Gpts+1:ihi,:,:), Gpts*lyTpts*NVAR_NUMBER, MPI_DOUBLE, &
                  r_blk, sendtag, comm, ierr)
    call MPI_RECV(solnVar(ihi+1:iGhi,:,:), Gpts*lyTpys*NVAR_NUMBER, MPI_DOUBLE, &
                  r_blk, recvtag, comm, status1, ierr)
  end if

end subroutine guardcell_fill
