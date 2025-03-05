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
      stop
    end subroutine
end module mpi_func
