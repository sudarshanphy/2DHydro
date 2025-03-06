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
  integer :: tag1 = 20, tag2 = 21
  integer :: xnum, ynum
  real, allocatable :: sendbuffer(:), receivebuffer(:)
  integer :: i, j, k, count

  solnVar(iGlo:, jGlo:, 1:) => mainVar(:,:,1:)

  xnum = Gpts*lyTpts*NVAR_NUMBER
  ynum = Gpts*lxTpts*NVAR_NUMBER

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
    sendbuffer(:) =  reshape(solnVar(ilo:ilo+Gpts-1,:,:), shape=[size(solnVar(ilo:ilo+Gpts-1,:,:))])

    call MPI_SENDRECV(sendbuffer, xnum, MPI_DOUBLE, l_blk, tag1, &
                   receivebuffer, xnum, MPI_DOUBLE, r_blk, tag1, &
                   comm, status1, ierr)
    if (r_blk /= MPI_PROC_NULL) solnVar(ihi+1:iGhi,:,:) = reshape(receivebuffer, shape=shape(solnVar(ihi+1:iGhi,:,:)))


    ! send data to right, receive from left
    sendbuffer(:) =  reshape(solnVar(ihi-Gpts+1:ihi,:,:), shape=[size(solnVar(ihi-Gpts+1:ihi,:,:))])
    call MPI_SENDRECV(sendbuffer, xnum, MPI_DOUBLE, r_blk, tag2, &
                   receivebuffer, xnum, MPI_DOUBLE, l_blk, tag2, &
                   comm, status1, ierr)
    if (l_blk /= MPI_PROC_NULL) solnVar(iGlo:ilo-1,:,:) = reshape(receivebuffer, shape=shape(solnVar(iGlo:ilo-1,:,:)))

    
    !print *, solnVar(iGlo:ilo-1,:,:)
    deallocate(sendbuffer, receivebuffer)
  endif
  if (yblk > 1) then
    ! communication in Y-dirextion
    allocate(sendbuffer(ynum))
    allocate(receivebuffer(ynum))
    
    ! send data to bottom, receive from top
    sendbuffer(:) =  reshape(solnVar(:,jlo:jlo+Gpts-1,:), shape=[size(solnVar(:,jlo:jlo+Gpts-1,:))])
    !count  = 1
    !do k = 1, NVAR_NUMBER
    !  do j = jlo, jlo+Gpts-1
    !    do i = iGlo, iGhi
    !      sendbuffer(count) = solnVar(i,j,k)
    !      count = count + 1
    !    end do
    !  end do
    !end do


    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, b_blk, tag2, &
                   receivebuffer, ynum, MPI_DOUBLE, t_blk, tag2, &
                   comm, status1, ierr)
    if (t_blk /= MPI_PROC_NULL) solnVar(:,jhi+1:jGhi,:) = reshape(receivebuffer, shape=shape(solnVar(:,jhi+1:jGhi,:)))
    !if (t_blk /= MPI_PROC_NULL) then
    !  count  = 1
    !  do k = 1, NVAR_NUMBER
    !    do j = jhi+1, jGhi
    !      do i = iGlo, iGhi
    !        solnVar(i,j,k) = receivebuffer(count)
    !        count = count + 1
    !      end do
    !    end do
    !  end do
    !end if


    ! send data to top, receive from bottom
    sendbuffer(:) =  reshape(solnVar(:,jhi-Gpts+1:jhi,:), shape=[size(solnVar(:,jhi-Gpts+1:jhi,:))])
    !count  = 1
    !do k = 1, NVAR_NUMBER
    !  do j = jhi-Gpts+1, jhi
    !    do i = iGlo, iGhi
    !      sendbuffer(count) = solnVar(i,j,k)
    !      count = count + 1
    !    end do
    !  end do
    !end do

    call MPI_SENDRECV(sendbuffer, ynum, MPI_DOUBLE, t_blk, tag1, &
                   receivebuffer, ynum, MPI_DOUBLE, b_blk, tag1, &
                   comm, status1, ierr)
    if (b_blk /= MPI_PROC_NULL) solnVar(:,jGlo:jlo-1,:) = reshape(receivebuffer, shape=shape(solnVar(:,jGlo:jlo-1,:)))
    !if (b_blk /= MPI_PROC_NULL) then
    !  count  = 1
    !  do k = 1, NVAR_NUMBER
    !    do j = jGlo, jlo-1
    !      do i = iGlo, iGhi
    !        solnVar(i,j,k) = receivebuffer(count)
    !        count = count + 1
    !      end do
    !    end do
    !  end do
    !end if

    
    deallocate(sendbuffer, receivebuffer)
  endif
  nullify(solnVar) 
end subroutine guardcell_fill
