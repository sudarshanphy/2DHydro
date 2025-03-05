#compiler
#FC=gfortran
FC=mpif90

#flags for the compiler
FFLAGS=-g -cpp -fdefault-real-8 

# Define a variable to indicate whether we are building for MHD
MHD ?= 0

# Library to include based on the target
#path for header file
ifeq ($(MHD), 1)
    HPATH=-I./src/MHD_h
else
    HPATH=-I./src/HD_h
endif

#files needed for compilation

#SRC=misc.f90 recon.f90 eos.f90 hllc.f90 test_hllc.f90
#SRC = io.f90 test_io.f90
SRC = sim_data.f90 mpi_func.f90 misc.f90 \
			read_par.f90 grid_func.f90 \
			eos.f90 glm.f90 \
			applyBC.f90 sim_init.f90 sim_restart.f90 io.f90 \
			recon.f90 get_flux.f90 get_resistive_flux.f90 rk2.f90 main.f90

#object files have .o extension
OBJ=${SRC:.f90=.o}

#Loop (%) over *.f90 files to create *.o files
#$@ represents the target files and $< represents the prerequisite files

%.o: src/%.f90
	$(FC) $(FFLAGS) $(HPATH) -c $< 

#executable to be created
EXEC=run_hd.exe
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) $(HPATH) $(OBJ) -o $@
	rm -f *.mod *.o  

#executable to be created
EXEC2=run_mhd.exe
$(EXEC2): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@ 
	rm -f *.mod *.o  

#clean all the compiled file and executable
clean:
	rm -f *.mod *.o
cleanall:
	rm -f *.mod *.o *.exe
cleandata:	
	rm ./output/*.dat

# Make executable for MHD run
mhd:
	$(MAKE) MHD=1 $(EXEC2)
