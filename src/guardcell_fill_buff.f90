subroutine guardcell_fill()
#include "param.h"
  use sim_data, only: ilo, jlo, ihi, jhi, &
                      iGlo, jGlo, iGhi, jGhi, &
                      mainVar, myrank, xblk, yblk, &
                      lxTpts, lyTpts, comm, ierr, status1, &
                      at_xlboundary, at_xrboundary, at_ylboundary, &
                      at_yrboundary, mainVar, Gpts

  use mpi
  implicit none
  real, pointer :: solnVar(:,:,:)
  integer :: l_blk, r_blk, t_blk, b_blk
  integer :: tag1 = 20, tag2 = 21, tag3 = 23, tag4 = 24
  integer :: xnum, ynum
  real, allocatable :: sendbufferl(:), sendbufferr(:), &
                    receivebufferl(:), receivebufferr(:)
  integer :: i, j, k, count, request

  solnVar(1:, iGlo:, jGlo:) => mainVar(:,:,:)

  xnum = Gpts*(jGhi - jGlo + 1)*NVAR_NUMBER
  ynum = Gpts*(iGhi - iGlo + 1)*NVAR_NUMBER

  ! see if we don't have neighbouring blocks
  l_blk = merge(MPI_PROC_NULL, myrank - 1, at_xlboundary)
  r_blk = merge(MPI_PROC_NULL, myrank + 1, at_xrboundary)
  t_blk = merge(MPI_PROC_NULL, myrank + xblk, at_yrboundary)
  b_blk = merge(MPI_PROC_NULL, myrank - xblk, at_ylboundary)

  !print *, "rank, l, r, t, b = ", myrank, l_blk, r_blk, t_blk, b_blk  

  if (xblk > 1) then
    ! communication in X-dirextion
    allocate(sendbufferl(xnum))
    allocate(receivebufferr(xnum))
    allocate(sendbufferr(xnum))
    allocate(receivebufferl(xnum))
    
    ! send data to left, receive from right 
    sendbufferl(:) =  reshape(solnVar(1:,ilo:ilo+Gpts-1,jGlo:jGhi), shape=[size(solnVar(1:,ilo:ilo+Gpts-1,jGlo:jGhi))])
    !if (l_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send left = ", myrank, solnVar(ilo:ilo+Gpts-1,jlo:jhi,DENS_VAR)
    !end if
    
    call MPI_SENDRECV(sendbufferl, xnum, MPI_DOUBLE, l_blk, tag1, &
                   receivebufferr, xnum, MPI_DOUBLE, r_blk, tag1, &
                   comm, status1, ierr)

    if (r_blk /= MPI_PROC_NULL) then
      solnVar(:,ihi+1:iGhi,jGlo:jGhi) = reshape(receivebufferr, shape=shape(solnVar(:,ihi+1:iGhi,jGlo:jGhi)))
      !print *, "myrank, recv right = ", myrank, solnVar(ihi+1:iGhi,jlo:jhi,DENS_VAR)
    end if

    !deallocate(sendbuffer, receivebuffer)
    ! send data to right, receive from left
    sendbufferr(:) =  reshape(solnVar(:,ihi-Gpts+1:ihi,jGlo:jGhi), shape=[size(solnVar(:,ihi-Gpts+1:ihi,jGlo:jGhi))])
    !if (r_blk /= MPI_PROC_NULL) then
    !  print *, "myank, send right = ", myrank, solnVar(ihi-Gpts+1:ihi,jlo:jhi,DENS_VAR)
    !end if
    call MPI_SENDRECV(sendbufferr, xnum, MPI_DOUBLE, r_blk, tag2, &
                   receivebufferl, xnum, MPI_DOUBLE, l_blk, tag2, &
                   comm, status1, ierr)

    if (l_blk /= MPI_PROC_NULL) then
      solnVar(:,iGlo:ilo-1,jGlo:jGhi) = reshape(receivebufferl, shape=shape(solnVar(:,iGlo:ilo-1,jGlo:jGhi)))
      !print *, "myrank, recv left = ", myrank, solnVar(iGlo:ilo-1,jlo:jhi, DENS_VAR)
    end if
    
    deallocate(sendbufferl, sendbufferr, receivebufferl, receivebufferr)
  endif
#if 0
  if (yblk > 1) then
    ! communication in Y-direction
    allocate(sendbuffer(ynum))
    allocate(receivebuffer(ynum))
    
    ! send data to bottom, receive from top
    sendbuffer(:) =  reshape(solnVar(:,iGlo:iGhi,jlo:jlo+Gpts-1), shape=[size(solnVar(:,iGlo:iGhi,jlo:jlo+Gpts-1))])
    !if (b_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send bottom = ", myrank, solnVar(ilo:ihi, jlo:jlo+Gpts-1, DENS_VAR)
    !end if

    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, b_blk, tag3, &
                   receivebuffer, ynum, MPI_DOUBLE, t_blk, tag3, &
                   comm, status1, ierr)
    if (t_blk /= MPI_PROC_NULL) then
      solnVar(:,iGlo:iGhi,jhi+1:jGhi) = reshape(receivebuffer, shape=shape(solnVar(:,iGlo:iGhi,jhi+1:jGhi)))
      !print *, "myrank, recv top = ", myrank, solnVar(ilo:ihi, jhi+1:jGhi, DENS_VAR)
    end if
    
    ! send data to top, receive from bottom
    sendbuffer(:) =  reshape(solnVar(1:,iGlo:iGhi,jhi-Gpts+1:jhi), shape=[size(solnVar(1:,iGlo:iGhi,jhi-Gpts+1:jhi))])
    !if (t_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send top = ", myrank, solnVar(ilo:ihi, jhi-Gpts+1:jhi, DENS_VAR)
    !end if

    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, t_blk, tag4, &
                   receivebuffer, ynum, MPI_DOUBLE, b_blk, tag4, &
                   comm, status1, ierr)
    if (b_blk /= MPI_PROC_NULL) then
      solnVar(:,iGlo:iGhi,jGlo:jlo-1) = reshape(receivebuffer, shape=shape(solnVar(:,iGlo:iGhi,jGlo:jlo-1)))
      !print *, "myrank, recv bottom = ", myrank, solnVar(ilo:ihi, jGlo:jlo-1,DENS_VAR)
    endif
    deallocate(sendbuffer, receivebuffer)
  endif
#endif
  nullify(solnVar) 
end subroutine guardcell_fill
