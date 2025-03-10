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
  real, allocatable :: sendbuffer(:), receivebuffer(:)
  integer :: i, j, k, count

  solnVar(iGlo:, jGlo:, 1:) => mainVar(:,:,1:)

  xnum = Gpts*(jGhi - jGlo + 1)*NVAR_NUMBER
  ynum = Gpts*(iGhi - iGlo + 1)*NVAR_NUMBER

  ! see if we don't have neighbouring blocks
  l_blk = merge(MPI_PROC_NULL, myrank - 1, at_xlboundary)
  r_blk = merge(MPI_PROC_NULL, myrank + 1, at_xrboundary)
  t_blk = merge(MPI_PROC_NULL, myrank + xblk, at_yrboundary)
  b_blk = merge(MPI_PROC_NULL, myrank - xblk, at_ylboundary)
 
  if (xblk > 1) then
    ! communication in X-dirextion
    allocate(sendbuffer(xnum))
    allocate(receivebuffer(xnum))
    
    ! send data to left, receive from right 
    sendbuffer(:) =  reshape(solnVar(ilo:ilo+Gpts-1,jGlo:jGhi,:), shape=[size(solnVar(ilo:ilo+Gpts-1,jGlo:jGhi,:))])
    !if (l_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send left = ", myrank, solnVar(ilo:ilo+Gpts-1,jlo:jhi,DENS_VAR)
    !end if
    
    call MPI_SENDRECV(sendbuffer, xnum, MPI_DOUBLE, l_blk, tag1, &
                   receivebuffer, xnum, MPI_DOUBLE, r_blk, tag1, &
                   comm, status1, ierr)
    if (r_blk /= MPI_PROC_NULL) then
      solnVar(ihi+1:iGhi,jGlo:jGhi,:) = reshape(receivebuffer, shape=shape(solnVar(ihi+1:iGhi,jGlo:jGhi,:)))
      !print *, "myrank, recv right = ", myrank, solnVar(ihi+1:iGhi,jlo:jhi,DENS_VAR)
    end if

    ! send data to right, receive from left
    sendbuffer(:) =  reshape(solnVar(ihi-Gpts+1:ihi,jGlo:jGhi,:), shape=[size(solnVar(ihi-Gpts+1:ihi,jGlo:jGhi,:))])
    !if (r_blk /= MPI_PROC_NULL) then
    !  print *, "myank, send right = ", myrank, solnVar(ihi-Gpts+1:ihi,jlo:jhi,DENS_VAR)
    !end if
    call MPI_SENDRECV(sendbuffer, xnum, MPI_DOUBLE, r_blk, tag2, &
                   receivebuffer, xnum, MPI_DOUBLE, l_blk, tag2, &
                   comm, status1, ierr)
    if (l_blk /= MPI_PROC_NULL) then
      solnVar(iGlo:ilo-1,jGlo:jGhi,:) = reshape(receivebuffer, shape=shape(solnVar(iGlo:ilo-1,jGlo:jGhi,:)))
      !print *, "myrank, recv left = ", myrank, solnVar(iGlo:ilo-1,jlo:jhi, DENS_VAR)
    end if
    
    !print *, solnVar(iGlo:ilo-1,:,:)
    deallocate(sendbuffer, receivebuffer)
  endif
  if (yblk > 1) then
    ! communication in Y-direction
    allocate(sendbuffer(ynum))
    allocate(receivebuffer(ynum))
    
    ! send data to bottom, receive from top
    sendbuffer(:) =  reshape(solnVar(iGlo:iGhi,jlo:jlo+Gpts-1,1:), shape=[size(solnVar(iGlo:iGhi,jlo:jlo+Gpts-1,1:))])
    !if (b_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send bottom = ", myrank, solnVar(ilo:ihi, jlo:jlo+Gpts-1, DENS_VAR)
    !end if

    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, b_blk, tag3, &
                   receivebuffer, ynum, MPI_DOUBLE, t_blk, tag3, &
                   comm, status1, ierr)
    if (t_blk /= MPI_PROC_NULL) then
      solnVar(iGlo:iGhi,jhi+1:jGhi,1:) = reshape(receivebuffer, shape=shape(solnVar(iGlo:iGhi,jhi+1:jGhi,1:)))
      !print *, "myrank, recv top = ", myrank, solnVar(ilo:ihi, jhi+1:jGhi, DENS_VAR)
    end if
    
    ! send data to top, receive from bottom
    sendbuffer(:) =  reshape(solnVar(iGlo:iGhi,jhi-Gpts+1:jhi,1:), shape=[size(solnVar(iGlo:iGhi,jhi-Gpts+1:jhi,1:))])
    !if (t_blk /= MPI_PROC_NULL) then
    !  print *, "myrank, send top = ", myrank, solnVar(ilo:ihi, jhi-Gpts+1:jhi, DENS_VAR)
    !end if

    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, t_blk, tag4, &
                   receivebuffer, ynum, MPI_DOUBLE, b_blk, tag4, &
                   comm, status1, ierr)
    if (b_blk /= MPI_PROC_NULL) then
      solnVar(iGlo:iGhi,jGlo:jlo-1,1:) = reshape(receivebuffer, shape=shape(solnVar(iGlo:iGhi,jGlo:jlo-1,:)))
      !print *, "myrank, recv bottom = ", myrank, solnVar(ilo:ihi, jGlo:jlo-1,DENS_VAR)
    endif
    deallocate(sendbuffer, receivebuffer)
  endif
  nullify(solnVar) 
end subroutine guardcell_fill
