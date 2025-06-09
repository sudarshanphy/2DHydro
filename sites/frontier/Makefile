#compiler
FC=ftn

#flag for main runs
FFLAGS=-fopenmp -s real64 -s integer32 -eI -eZ -ef -f free 

#profiling
#rocprof --stats --hsa-trace --hip-trace ./exec.exe

#files needed for compilation

SRC = sim_data.f90 mpi_func.f90 misc.f90 \
			read_par.f90 grid_func.f90 guardcell_fill.f90\
			eos.f90 glm.f90 \
			applyBC.f90 sim_init.f90 sim_restart.f90 io.f90 \
			recon.f90 get_flux.f90 get_resistive_flux.f90 rk2.f90 main.f90

#object files have .o extension
OBJ=${SRC:.f90=.o}

#Loop (%) over *.f90 files to create *.o files
#$@ represents the target files and $< represents the prerequisite files

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

#executable to be created
EXEC=run.exe
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@
	@echo "Command to run on Frontier:"
	@echo "stdbuf -o0 -e0 srun -N Nnodes -n Nprocs --ntasks-per-gpu=Ngpus ./$(EXEC)"

#clean all the compiled file and executable
clean:
	rm -f *.mod *.o *.i
