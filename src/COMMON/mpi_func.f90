module mpi_func
#include "param.h"
  implicit none
  contains
  
    subroutine init_procs()
      use sim_data, only: comm, myrank, nprocs, &
                          ierr
      use mpi 
      implicit none

      comm = MPI_COMM_WORLD ! used frequently
      ! initialize mpi
      call MPI_INIT(ierr)
      
      ! get number of processors
      call MPI_COMM_SIZE(comm, nprocs, ierr)
      
      ! get processor rank
      call MPI_COMM_RANK(comm, myrank, ierr)

    end subroutine
    
    subroutine check_procs()
      use sim_data, only: comm, myrank, nprocs, &
                          ierr, xblk, yblk, errcode
      use mpi 
      implicit none

      if (nprocs /= xblk * yblk) then
        print *, "nprocs, xblk, yblk = ", nprocs, xblk, yblk
         if (myrank == MASTER_PROC) then
           print *, "number of processors should be = xblk * yblk"
         end if
         call MPI_Abort(comm, errcode, ierr)
         stop
      end if

    end subroutine

    subroutine finalize_procs()
      use sim_data, only: ierr
      use mpi
      implicit none
      ! finalize mpi
      call MPI_FINALIZE(ierr)
      !stop   !print IEEE meesages at the botom if not commented
    end subroutine

    subroutine create_mpiDatatypes()
      use sim_data, only: Gpts, ierr, dtype_mpi, &
                          iGlo, iGhi, jGlo, jGhi
      use mpi
      implicit none
      integer, dimension(1:NDIM) :: blk_extent
      integer :: exch1, exch2

      blk_extent(1) = iGhi - iGlo + 1
      blk_extent(2) = jGhi - jGlo + 1

      ! for top and bottom, ydir
      call MPI_TYPE_CONTIGUOUS(blk_extent(1)*NVAR_NUMBER*Gpts, &
                               MPI_DOUBLE, exch2, ierr)
      call MPI_TYPE_COMMIT(exch2,ierr)

      ! for left and right, xdir
      CALL MPI_TYPE_VECTOR(blk_extent(2), NVAR_NUMBER*Gpts, &
               NVAR_NUMBER*blk_extent(1), MPI_DOUBLE, exch1, ierr)
      CALL MPI_TYPE_COMMIT(exch1, ierr)

      dtype_mpi(1) = exch1
      dtype_mpi(2) = exch2

    end subroutine
end module mpi_func
